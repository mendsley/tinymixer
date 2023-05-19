/*-
 * Copyright 2011-2013 Matthew Endsley
 * All rights reserved
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted providing that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef TINYMIXER_MIXER_H
#define TINYMIXER_MIXER_H

#include <stdint.h>

struct tinymixer_buffer;
struct tinymixer_channel;

struct tinymixer_callbacks
{
	// user context. tinymixer will not modify or interpret this value. It is
	// simply passed as-is to the callback functions
	void* opaque = nullptr;

	// allocate and free memory. If null, tinymixer will use malloc/free
	void* (*allocate)(void* opaque, int bytes) = nullptr;
	void  (*free)(void* opaque, void* pointer) = nullptr;

	// allow the client to add additional audio into the mix before applying
	// the effects pass
	void (*pre_effects)(void* opaque, float* samples, int nsamples, float gain) = nullptr;

	// called when a channel stops playing
	void (*channel_complete)(void* opaque, void* channel_opaque, tinymixer_channel channel) = nullptr;
};

struct tinymixer_buffer_callbacks
{
	void (*on_destroy)(void* buffer_opqaue);

	void* (*start_source)(void* buffer_opaque);
	void (*end_source)(void* buffer_opaque, void* source_opaque);
	int (*request_samples)(void* buffer_opaque, void* source_opaque, const float** left, const float** right, int nsamples);
};

void tinymixer_init(tinymixer_callbacks callbacks, int sample_rate);
void tinymixer_getsamples(float* samples, int nsamples);
void tinymixer_set_mastergain(float gain);

void tinymixer_create_buffer_interleaved_s16le(int channels, const int16_t* pcm_data, int pcm_data_size
		, const tinymixer_buffer** handle);
void tinymixer_create_buffer_interleaved_float(int channels, const float* pcm_data, int pcm_data_size
		, const tinymixer_buffer** handle);
void tinymixer_create_buffer_vorbis_stream(const void* data, int ndata
		, void* opaque, void (*closed)(void*), const tinymixer_buffer** handle);
void tinymixer_create_buffer_custom_stream(void* opaque, tinymixer_buffer_callbacks callbacks, const tinymixer_buffer** buffer);
int tinymixer_get_buffer_size(const tinymixer_buffer* handle);
void tinymixer_release_buffer(const tinymixer_buffer* handle);

void tinymixer_update_listener(const float* position);
void tinymixer_set_base_gain(int index, float gain);
void tinymixer_set_callback_gain(float gain);
void tinymixer_effects_compressor(const float thresholds[2], const float multipliers[2]
		, float attack_seconds, float release_seconds);

bool tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch
		, tinymixer_channel* channel);
bool tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, const float* position
		, float distance_min, float distance_max, tinymixer_channel* channel);

bool tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch
		, tinymixer_channel* channel);
bool tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch
		, const float* position, float distance_min, float distance_max, tinymixer_channel* channel);

void tinymixer_channel_set_opaque(tinymixer_channel channel, void* opaque);
void tinymixer_channel_stop(tinymixer_channel channel);
void tinymixer_channel_set_position(tinymixer_channel channel, const float* position);
void tinymixer_channel_fadeout(tinymixer_channel channel, float seconds);
void tinymixer_channel_set_gain(tinymixer_channel channel, float gain);
void tinymixer_channel_set_frequency(tinymixer_channel channel, float frequency);

float tinymixer_channel_get_gain(tinymixer_channel channel);
void tinymixer_stop_all_sources();

struct tinymixer_channel
{
	int index = 0;
};

static inline bool tinymixer_channel_isvalid(tinymixer_channel channel) { return channel.index != 0; }
static inline bool operator==(tinymixer_channel lhs, tinymixer_channel rhs) { return lhs.index == rhs.index; }

struct tinymixer_resampler
{
	float ideal_rate;
	float prev_samples[2];
};

void tinymixer_resampler_init(tinymixer_resampler* resampler, int input_sample_rate, int output_sample_rate);
void tinymixer_resampler_init_rate(tinymixer_resampler* resampler, float ideal_rate);
int tinymixer_resampler_calculate_input_samples(const tinymixer_resampler* resampler, int output_samples);
int tinymixer_resampler_calculate_output_samples(const tinymixer_resampler* resampler, int input_samples);
void tinymixer_resample_stereo(tinymixer_resampler* resampler, const float* input, int num_input_samples, float* output, int num_output_samples);
void tinymixer_resample_mono(tinymixer_resampler* resampler, const float* input, int num_input_samples, float* output, int num_output_samples);

struct tinymixer_lowpass_filter
{
	float cutoff_frequency; // 1250.0 is a decent starting point
	float sample_rate;

	float channel_history[2];
};

void tinymixer_lowpass_filter_init(tinymixer_lowpass_filter* filter, float cutoff_frequency, float sample_rate);
void tinymixer_lowpass_filter_apply(tinymixer_lowpass_filter* filter, float* output, float* input, int num_samples, int num_channels);

#endif
