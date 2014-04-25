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
typedef uint8_t tinymixer_loop;
typedef void (*tinymixer_callback)(float* samples, int nsamples, float gain);

void tinymixer_init(int sample_rate, tinymixer_callback callback);
void tinymixer_getsamples(float* samples, int nsamples);
void tinymixer_set_mastergain(float gain);

void tinymixer_create_buffer_interleaved_s16le(int channels, const int16_t* pcm_data, int pcm_data_size, const tinymixer_buffer** handle);
void tinymixer_create_buffer_interleaved_float(int channels, const float* pcm_data, int pcm_data_size, const tinymixer_buffer** handle);
tinymixer_buffer* tinymixer_create_buffer(int channels, int nsamples);
void tinymixer_buffer_getchannelbuffers(tinymixer_buffer* handle, float** channels);
int tinymixer_get_buffer_size(const tinymixer_buffer* handle);
void tinymixer_release_buffer(const tinymixer_buffer* handle);

void tinymixer_update_listener(const float* position);
void tinymixer_set_base_gain(int index, float gain);
void tinymixer_set_callback_gain(float gain);
void tinymixer_effects_compressor(const float thresholds[2], const float multipliers[2], float attack_seconds, float release_seconds);

void tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch);
void tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, const float* position, float distance_min, float distance_max);

void tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, tinymixer_loop* loop);
void tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch,
		const float* position, float distance_min, float distance_max, tinymixer_loop* loop);
void tinymixer_remove_loop(tinymixer_loop loop);
void tinymixer_loop_set_position(tinymixer_loop loop, const float* position);
void tinymixer_loop_fadeout(tinymixer_loop loop, float seconds);
void tinymixer_loop_set_gain(tinymixer_loop loop, float gain);
void tinymixer_loop_set_frequency(tinymixer_loop loop, float frequency);

float tinymixer_loop_get_gain(tinymixer_loop loop);
void tinymixer_stop_all_sources();

#endif
