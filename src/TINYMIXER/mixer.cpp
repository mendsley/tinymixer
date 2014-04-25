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

#include <TINYMIXER/mixer.h>
#include <stdint.h>

#ifndef tinymixer_sqrtf
#include <math.h>
#define tinymixer_sqrtf sqrtf
#endif

#ifndef tinymixer_abs
#include <stdlib.h>
#define tinymixer_abs abs
#endif

#ifndef tinymixer_fabs
#include <math.h>
#define tinymixer_fabs fabs
#endif

#if !defined(tinymixer_free) || !defined(tinymixer_alloc)
#include <stdlib.h>
#define tinymixer_alloc malloc
#define tinymixer_free free
#endif

#if !defined(tinymixer_memset)
#include <string.h>
#define tinymixer_memset memset
#endif

namespace {
struct SourceFlags {
	enum Enum {
		Playing		= 1 << 0,
		Positional	= 1 << 1,
		Looping		= 1 << 2,
		FadeOut		= 1 << 3,
		Frequency	= 1 << 4,
	};
};

struct Buffer;
struct Source;

struct buffer_functions {
	void (*on_destroy)(Buffer* buffer);

	void (*start_source)(Source* source);
	void (*end_source)(Source* source);
	int (*request_samples)(Source* source, const float** left, const float** right, int nsamples);

	int (*get_buffer_size)(const Buffer* buffer);
};

struct Buffer {
	const buffer_functions *funcs;
	int32_t refcnt;
	int32_t nsamples;
	uint8_t nchannels;
	// float smaples[nsamples*nschannels];
};

struct static_source_data {
	int32_t sample_pos;
};

struct Source {
	const Buffer* buffer;
	union {
		static_source_data static_source;
	} instance_data;
	float position[3];
	float fadeout_per_sample;
	float gain_base;
	float distance_min;
	float distance_difference;
	float frequency;
	uint16_t frame_age;
	uint8_t flags;
	uint8_t gain_base_index;
};

static const int c_ngaintypes = 8;
static const int c_nsources = 32;
static const int c_nsamples = 2048;
static const float c_fnsamples = (float)c_nsamples;
static float c_speakerdist = 0.17677669529663688110021109052621f; // 1/(4 *sqrtf(2))

struct Mixer {
	tinymixer_callback callback;
	float position[3];
	float gain_master;
	float gain_base[c_ngaintypes];
	float gain_callback;
	int32_t sample_rate;
	float compressor_last_samples[2];
	float compressor_thresholds[2];
	float compressor_multipliers[2];
	float compressor_factor;
	float compressor_attack_per1ksamples;
	float compressor_release_per1ksamples;
	int32_t samples_remaining;
	Source sources[c_nsources];
	float buffer[2*c_nsamples];
	float scratch[2*c_nsamples];
};
}

static const int c_quantize = 1024;
static const int c_quantize_mask = (c_quantize - 1);
static Mixer g_mixer;

static inline float mixer_clamp(float v, float min, float max) { return min > v ? min : (v > max ? max : v); }
static inline int32_t mixer_clamp(int32_t v, int32_t min, int32_t max) { return min > v ? min : (v > max ? max : v); }
static inline int mixer_min(int a, int b) { return (a < b) ? a : b; }
static inline int mixer_max(int a, int b) { return (b < a) ? a : b; }
static inline float mixer_dist(const float* a, const float* b) {
	const float distsq = (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + (a[2] - b[2])*(a[2] - b[2]);
	return tinymixer_sqrtf(distsq);
}
static inline void mixer_vcopy(float* v, const float* a) { v[0] = a[0], v[1] = a[1], v[2] = a[2]; }

static void addref(Buffer* buffer) {
	++buffer->refcnt;
}

static void decref(Buffer* buffer) {
	if (0 == --buffer->refcnt) {
		buffer->funcs->on_destroy(buffer);
		tinymixer_free(buffer);
	}
}

static void kill_source(Source* source) {
	source->buffer->funcs->end_source(source);

	decref((Buffer*)source->buffer);
	source->buffer = 0;
	source->flags = 0;
}

static Source *find_source() {
	Source* best_source = 0;
	uint16_t best_age = 0xffff;

	for (int ii = 0; ii < c_nsources; ++ii) {
		Source* source = &g_mixer.sources[ii];
		if (!source->buffer)
			return source;

		if (0 == (source->flags & SourceFlags::Looping)) {
			const uint16_t age = source->frame_age;
			if (age < best_age) {
				best_source = source;
				best_age = age;
			}
		}
	}

	if (NULL != best_source) {
		if (best_source->buffer) {
			kill_source(best_source);
		}
	}
	return best_source;
}

static void resample_mono(float* out, int nout, const float* in, int qfreq) {
	for (int qpos = 0; nout; --nout, qpos += qfreq) {
		const int qindex = qpos / c_quantize;
		const int qinterp = qpos % c_quantize;

		*out++ = in[qindex] + (float)qinterp * (in[qindex + 1] - in[qindex]) / (float)c_quantize;
	}
}

static void render(Source* source, float* buffer, const float gain[2]) {

	float* left = buffer;
	float* right = buffer + c_nsamples;

	int remaining = c_nsamples;
	while (remaining) {
		int samples_read = remaining;
		int samples_written = samples_read;

		const float* srcleft;
		const float* srcright;

		// source has a non-1.0f frequency shift
		if (source->flags & SourceFlags::Frequency) {
			const int qfreq = (int)(source->frequency * c_quantize);
			samples_read = source->buffer->funcs->request_samples(source, &srcleft, &srcright, mixer_min((samples_written * qfreq) / c_quantize, c_nsamples));
			if (samples_read == 0) {
				// source is no longer playing
				source->flags &= ~SourceFlags::Playing;
				return;
			}

			samples_written = (samples_read * c_quantize) / qfreq;
			if (samples_written == 0)
				return;

			// resample source into scratch space
			resample_mono(g_mixer.scratch, samples_written, srcleft, qfreq);
			if (srcleft == srcright) {
				srcleft = srcright = g_mixer.scratch;
			} else {
				// resample stereo channel
				resample_mono(g_mixer.scratch + samples_written, samples_written, srcright, qfreq);
				srcleft = g_mixer.scratch;
				srcright = g_mixer.scratch + samples_written;
			}

		} else {
			samples_read = source->buffer->funcs->request_samples(source, &srcleft, &srcright, samples_read);
			if (samples_read == 0) {
				// source is no longer playing
				source->flags &= ~SourceFlags::Playing;
				return;
			}
			samples_written = samples_read;
		}

		// render the source to the output mix
		for (int ii = 0; ii < samples_written; ++ii)
			*left++ += gain[0] * srcleft[ii];
		for (int ii = 0; ii < samples_written; ++ii)
			*right++ += gain[1] * srcright[ii];

		remaining -= samples_written;
	}
}

static void render_effects(float* buffer) {
	float compressor_factor = g_mixer.compressor_factor;

	// get maximum absolute power level from the rendered buffer, and adjust the compressor factor
	float max_power = 0;
	for (int ii = 0; ii < c_nsamples; ++ii) {
		const float power = tinymixer_fabs(buffer[ii]);
		if (power > max_power)
			max_power = power;
	}

	float target_compressor_factor = 1.0f;
	if (max_power > g_mixer.compressor_thresholds[1])
		target_compressor_factor = g_mixer.compressor_multipliers[1];
	else if (max_power > g_mixer.compressor_thresholds[0])
		target_compressor_factor = g_mixer.compressor_multipliers[0];

	float attack_release = 1.0f;
	if (target_compressor_factor < compressor_factor)
		attack_release = g_mixer.compressor_attack_per1ksamples;
	else if (target_compressor_factor > compressor_factor)
		attack_release = g_mixer.compressor_release_per1ksamples;

	// linearly interp compressor_factor toward the target compressor value
	const float interp = attack_release * c_fnsamples;
	compressor_factor = compressor_factor + interp*(target_compressor_factor - compressor_factor);
	compressor_factor = mixer_clamp(compressor_factor, g_mixer.compressor_multipliers[1], 1.0f);

    // 2-pass compressor to limit dynamic range of audio clipping levels
	if (compressor_factor < 1.0f) {
		for (int cc = 0; cc < 2; ++cc) {
			float prev_sample = g_mixer.compressor_last_samples[cc];

			float sample = 0;
			float* channel = buffer + cc*c_nsamples;
			for (int ii = 0; ii < c_nsamples; ++ii) {
				float sample = channel[ii];

				// uhhh... linear space? really??
				float diff = sample - prev_sample;
				sample = prev_sample + compressor_factor*diff;
				channel[ii] = sample;
			}

			g_mixer.compressor_last_samples[cc] = sample;
		}
	}

	g_mixer.compressor_factor = compressor_factor;
}

static void mix(float* buffer) {
	int nplaying = 0;
	int playing[c_nsources];
	float gain[c_nsources][2];

	// build active sounds
	for (int ii = 0; ii < c_nsources; ++ii) {
		const Source* source = &g_mixer.sources[ii];
		if (source->flags & SourceFlags::Playing) {
			playing[nplaying] = ii;
			++nplaying;
		}
	}

	// Update source gains
	for (int ii = 0; ii < nplaying; ++ii) {
		const Source* source = &g_mixer.sources[playing[ii]];

		const float gain_base = g_mixer.gain_master * g_mixer.gain_base[source->gain_base_index] * source->gain_base;
		gain[ii][0] = gain_base;
		gain[ii][1] = gain_base;

		if (source->flags & SourceFlags::Positional) {
			const float dist = mixer_dist(g_mixer.position, source->position);
			if (dist > 1.0e-8f) {
				// attenuation factor due to distance
				const float gain_distance = 1.0f - (dist - source->distance_min) / source->distance_difference;

				// attenuation factor due to panning (position audio)
				const float gain_panning = (source->position[0] - g_mixer.position[0]) / dist;

				gain[ii][0] *= gain_distance * (1.0f + gain_panning * -c_speakerdist);
				gain[ii][1] *= gain_distance * (1.0f + gain_panning * +c_speakerdist);
			}
		}

		// clamp gains
		gain[ii][0] = mixer_clamp(gain[ii][0], 0.0f, 1.0f);
		gain[ii][1] = mixer_clamp(gain[ii][1], 0.0f, 1.0f);
	}

	tinymixer_memset(buffer, 0, sizeof(float)*2*c_nsamples);

	// render playing sources
	for (int ii = 0; ii < nplaying; ++ii) {
		Source* source = &g_mixer.sources[playing[ii]];
		render(source, buffer, gain[ii]);
	}

	// allow application to apply a premixed track (such as music)
	if (g_mixer.callback)
		(g_mixer.callback)(buffer, c_nsamples, g_mixer.gain_callback);

	// render effects
	render_effects(buffer);

	// perform source-level post processing
	for (int ii = 0; ii < nplaying; ++ii) {
		Source* source = &g_mixer.sources[playing[ii]];
		++source->frame_age;

		// handle fadeout->stop
		if (source->flags & SourceFlags::FadeOut) {
			source->gain_base -= source->fadeout_per_sample * c_nsamples;
			if (source->gain_base <= 0.0f) {
				source->flags = 0;
			}
		}
	}

	// cleanup dead sources
	for (int ii = 0; ii < nplaying; ++ii) {
		Source* source = &g_mixer.sources[playing[ii]];

		if (0 == (source->flags & SourceFlags::Playing)) {
			kill_source(source);
		}
	}
}

void tinymixer_getsamples(float* samples, int nsamples) {
	// was data leftover after the previous call to getsamples? Copy that out here
	while (nsamples && g_mixer.samples_remaining) {
		const int samples_to_mix = mixer_min(nsamples, g_mixer.samples_remaining);
		const int offset = c_nsamples - g_mixer.samples_remaining;

		// clip and interleave
		for (int cc = 0; cc < 2; ++cc)
			for (int ii = 0; ii < samples_to_mix; ++ii)
				samples[cc + 2*ii] = mixer_clamp(g_mixer.buffer[cc*c_nsamples + offset + ii], -1.0f, 1.0f);

		g_mixer.samples_remaining -= samples_to_mix;
		samples += (2*samples_to_mix);
		nsamples -= samples_to_mix;
	}

	// Copy out samples
	if (nsamples) {
		const int samples_to_mix = mixer_min(nsamples, c_nsamples);
		mix(g_mixer.buffer);
		g_mixer.samples_remaining = c_nsamples;

		// clip and interleave
		for (int cc = 0; cc < 2; ++cc)
			for (int ii = 0; ii < samples_to_mix; ++ii)
				samples[cc + 2*ii] = mixer_clamp(g_mixer.buffer[cc*c_nsamples + ii], -1.0f, 1.0f);


		g_mixer.samples_remaining -= samples_to_mix;
		samples += (2*samples_to_mix);
		nsamples -= samples_to_mix;
	}
}

void tinymixer_set_mastergain(float gain) {
	g_mixer.gain_master = gain;
}

static Source* add(const tinymixer_buffer* handle, int gain_index, float gain) {
	Source* source = find_source();
	if (!source)
		return 0;

	source->buffer = (const Buffer*)handle;
	source->gain_base = gain;
	source->gain_base_index = (uint8_t)gain_index;
	source->frame_age = 0;

	source->buffer->funcs->start_source(source);

	addref((Buffer*)handle);
	return source;
}

static void play(Source* source) {
	source->flags |= SourceFlags::Playing;
}


static void static_sample_buffer_on_destroy(Buffer* buffer) {
}

static void static_sample_buffer_start_source(Source* source) {
	source->instance_data.static_source.sample_pos = 0;
}

static void static_sample_buffer_end_source(Source*source) {
}

static int static_sample_buffer_request_samples(Source* source, const float** left, const float** right, int nsamples) {
	int sample_pos = source->instance_data.static_source.sample_pos;

	// handle looping sources
	if (sample_pos == source->buffer->nsamples) {
		if (source->flags & SourceFlags::Looping) {
			sample_pos = 0;
		} else {
			return 0;
		}
	}

	nsamples = mixer_min(source->buffer->nsamples - sample_pos, nsamples);
	const float* srcleft = (float*)(source->buffer + 1) + sample_pos;

	*left = srcleft;
	if (source->buffer->nchannels == 1)
		*right = srcleft;
	else {
		*right = srcleft + source->buffer->nsamples;
	}

	source->instance_data.static_source.sample_pos += nsamples;
	return nsamples;
}

static int static_sample_buffer_get_buffer_size(const Buffer* buffer) {
	return sizeof(Buffer) + sizeof(float)*buffer->nchannels*buffer->nsamples;
}

static buffer_functions static_sample_functions = {
	static_sample_buffer_on_destroy,
	static_sample_buffer_start_source,
	static_sample_buffer_end_source,
	static_sample_buffer_request_samples,
	static_sample_buffer_get_buffer_size,
};

void tinymixer_create_buffer_interleaved_s16le(int channels, const int16_t* pcm_data, int pcm_data_size, const tinymixer_buffer** handle) {
	const int nsamples = pcm_data_size/sizeof(uint16_t)/channels;

	Buffer* buffer = (Buffer*)tinymixer_alloc(sizeof(Buffer) + nsamples*channels*sizeof(float));
	buffer->funcs = &static_sample_functions;
	buffer->refcnt = 1;
	buffer->nchannels = (uint8_t)channels;
	buffer->nsamples = nsamples;

	// copy samples
	const int16_t* source = (const int16_t*)pcm_data;
	float* dest = (float*)(buffer + 1);
	for (int cc = 0; cc < channels; ++cc) {
		for (int ii = 0; ii < nsamples; ++ii)
			*dest++ = (float)source[channels*ii + cc] / (float)0x8000;
	}

	*handle = (tinymixer_buffer*)buffer;
}

void tinymixer_create_buffer_interleaved_float(int channels, const float* pcm_data, int pcm_data_size, const tinymixer_buffer** handle) {
	const int nsamples = pcm_data_size/sizeof(float)/channels;

	Buffer* buffer = (Buffer*)tinymixer_alloc(sizeof(Buffer) + nsamples*channels*sizeof(float));
	buffer->funcs = &static_sample_functions;
	buffer->refcnt = 1;
	buffer->nchannels = (uint8_t)channels;
	buffer->nsamples = nsamples;

	// copy samples
	const float* source = (const float*)pcm_data;
	float* dest = (float*)(buffer + 1);
	for (int cc = 0; cc < channels; ++cc) {
		for (int ii = 0; ii < nsamples; ++ii)
			*dest++ = source[channels*ii + cc];//(int16_t)mixer_clamp((int32_t)(source[channels*ii + cc] * (float)0x8000), (int16_t)0x8000, 0x7fff);
	}

	*handle = (tinymixer_buffer*)buffer;
}

int tinymixer_get_buffer_size(const tinymixer_buffer* handle) {
	const Buffer* buffer = (const Buffer*)handle;
	return buffer->funcs->get_buffer_size(buffer);
}

void tinymixer_release_buffer(const tinymixer_buffer* handle) {
	decref((Buffer*)handle);
}

void tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch) {
	Source* source = add(handle, gain_index, gain);
	if (source) {
		if (pitch != 1.0f) {
			// clear frequency shift if ~0.0f
			const float diff = pitch - 1.0f;
			if (diff*diff < 1.0e-8f) {
				source->flags &= ~SourceFlags::Frequency;
			} else {
				source->frequency = pitch;
				source->flags |= SourceFlags::Frequency;
			}
		} else {
			source->flags &= ~SourceFlags::Frequency;
			source->frequency = 1.0f;
		}
		play(source);
	}
}

void tinymixer_add(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, const float* position, float distance_min, float distance_max) {
	Source *source = add(handle, gain_index, gain);
	if (source) {
		if (pitch != 1.0f) {
			// clear frequency shift if ~0.0f
			const float diff = pitch - 1.0f;
			if (diff*diff < 1.0e-8f) {
				source->flags &= ~SourceFlags::Frequency;
			} else {
				source->frequency = pitch;
				source->flags |= SourceFlags::Frequency;
			}
		} else {
			source->flags &= ~SourceFlags::Frequency;
			source->frequency = 1.0f;
		}
		mixer_vcopy(source->position, position);
		source->flags |= SourceFlags::Positional;
		source->distance_min = distance_min;
		source->distance_difference = (distance_max - distance_min);
		play(source);
	}
}

static Source* add_loop(const tinymixer_buffer* handle, int gain_index, float gain, tinymixer_loop* loop) {
	Source* source = add(handle, gain_index, gain);
	if (!source) {
		*loop = 0;
		return 0;
	}

	source->flags |= SourceFlags::Looping;
	*loop = (tinymixer_loop)(source - g_mixer.sources) + 1;
	return source;
}

void tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, tinymixer_loop* loop) {
	Source* source = add_loop(handle, gain_index, gain, loop);
	if (source)
	{
		if (pitch != 1.0f) {
			// clear frequency shift if ~0.0f
			const float diff = pitch - 1.0f;
			if (diff*diff < 1.0e-8f) {
				source->flags &= ~SourceFlags::Frequency;
			} else {
				source->frequency = pitch;
				source->flags |= SourceFlags::Frequency;
			}
		} else {
			source->flags &= ~SourceFlags::Frequency;
			source->frequency = 1.0f;
		}
		play(source);
	}
}

void tinymixer_add_loop(const tinymixer_buffer* handle, int gain_index, float gain, float pitch, const float* position, float distance_min, float distance_max, tinymixer_loop* loop) {
	Source* source = add_loop(handle, gain_index, gain, loop);
	if (source) {
		mixer_vcopy(source->position, position);
		source->flags |= SourceFlags::Positional;
		source->distance_min = distance_min;
		source->distance_difference = (distance_max - distance_min);
		if (pitch != 1.0f) {
			// clear frequency shift if ~0.0f
			const float diff = pitch - 1.0f;
			if (diff*diff < 1.0e-8f) {
				source->flags &= ~SourceFlags::Frequency;
			} else {
				source->frequency = pitch;
				source->flags |= SourceFlags::Frequency;
			}
		} else {
			source->flags &= ~SourceFlags::Frequency;
			source->frequency = 1.0f;
		}
		play(source);
	}
}

void tinymixer_remove_loop(tinymixer_loop loop) {
	Source* source = &g_mixer.sources[loop - 1];
	source->flags &= ~(SourceFlags::Playing | SourceFlags::Looping);
}

void tinymixer_loop_set_position(tinymixer_loop loop, const float* position) {
	Source* source = &g_mixer.sources[loop - 1];
	mixer_vcopy(source->position, position);
}

void tinymixer_loop_fadeout(tinymixer_loop loop, float seconds) {
	Source* source = &g_mixer.sources[loop - 1];
	source->fadeout_per_sample = 1.0f / (seconds * g_mixer.sample_rate);
	source->flags |= SourceFlags::FadeOut;
}

void tinymixer_loop_set_gain(tinymixer_loop loop, float gain) {
	Source* source = &g_mixer.sources[loop - 1];
	source->gain_base = gain;
	source->flags &= ~SourceFlags::FadeOut;
}

float tinymixer_loop_get_gain(tinymixer_loop loop) {
	Source* source = &g_mixer.sources[loop - 1];
	return source->gain_base;
}

void tinymixer_loop_set_frequency(tinymixer_loop loop, float frequency) {
	Source* source = &g_mixer.sources[loop - 1];

	// clear frequency shift if ~0.0f
	const float diff = frequency - 1.0f;
	if (diff*diff < 1.0e-8f) {
		source->flags &= ~SourceFlags::Frequency;
	} else {
		source->frequency = frequency;
		source->flags |= SourceFlags::Frequency;
	}
}

void tinymixer_init(int sample_rate, tinymixer_callback callback) {
	g_mixer.gain_master = 1.0f;
	for (int ii = 0; ii < c_ngaintypes; ++ii)
		g_mixer.gain_base[ii] = 1.0f;

	g_mixer.sample_rate = sample_rate;
	g_mixer.callback = callback;
	g_mixer.samples_remaining = 0;

	const float default_thresholds[2] = {1.0f, 1.0f};
	const float default_multipliers[2] = {1.0f, 1.0f};
	const float default_attack = 0.0f;
	const float default_release = 0.0f;
	tinymixer_effects_compressor(default_thresholds, default_multipliers, default_attack, default_release);
	g_mixer.compressor_factor = 1.0f;
	g_mixer.compressor_last_samples[0] = g_mixer.compressor_last_samples[1] = 0;
}

void tinymixer_update_listener(const float* position) {
	mixer_vcopy(g_mixer.position, position);
}

void tinymixer_set_base_gain(int index, float gain) {
	g_mixer.gain_base[index] = gain;
}

void tinymixer_set_callback_gain(float gain) {
	g_mixer.gain_callback = gain;
}

void tinymixer_effects_compressor(const float thresholds[2], const float multipliers[2], float attack_seconds, float release_seconds) {
	g_mixer.compressor_thresholds[0] = mixer_clamp(thresholds[0], 0.0f, 1.0f);
	g_mixer.compressor_thresholds[1] = mixer_clamp(thresholds[1], 0.0f, 1.0f);
	g_mixer.compressor_multipliers[0] = mixer_clamp(multipliers[0], 0.0f, 1.0f);
	g_mixer.compressor_multipliers[1] = mixer_clamp(multipliers[1], 0.0f, 1.0f);
	g_mixer.compressor_attack_per1ksamples = 1.0f / (attack_seconds * (float)g_mixer.sample_rate);
	g_mixer.compressor_release_per1ksamples = 1.0f / (release_seconds * (float)g_mixer.sample_rate);
}

void tinymixer_stop_all_sources() {
	for (int ii = 0; ii < c_nsources; ++ii) {
		Source* source = &g_mixer.sources[ii];
		if (source->buffer) {
			kill_source(source);
		}
	}
}
