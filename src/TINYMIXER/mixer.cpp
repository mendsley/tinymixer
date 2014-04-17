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

struct Buffer {
	int32_t refcnt;
	int32_t nsamples;
	uint8_t nchannels;
	// int16_t smaples[nsamples*nschannels];
};

struct Source {
	const Buffer* buffer;
	int32_t sample_pos;
	float position[3];
	float fadeout_per_sample;
	float gain_base;
	float distance_min;
	float distance_difference;
	float frequency;
	uint8_t flags;
	uint8_t gain_base_index;
};

static const int c_ngaintypes = 8;
static const int c_nsources = 32;
static const int c_nsamples = 2048;
static float c_speakerdist = 0.17677669529663688110021109052621f; // 1/(4 *sqrtf(2))

struct Mixer {
	tinymixer_callback callback;
	float position[3];
	float gain_master;
	float gain_base[c_ngaintypes];
	float gain_callback;
	int32_t sample_rate;
	int32_t compressor_last_samples[2];
	int32_t compressor_thresholds[2];
	int32_t compressor_multipliers[2];
	int32_t compressor_factor;
	int32_t compressor_attack_per1ksamples;
	int32_t compressor_release_per1ksamples;
	int32_t samples_remaining;
	Source sources[c_nsources];
	int32_t buffer[2*c_nsamples];
};
}

static const int c_quantize = 1024;
static const int c_quantize_mask = (c_quantize - 1);
static Mixer g_mixer;

static inline float mixer_clamp(float v, float min, float max) { return min > v ? min : (v > max ? max : v); }
static inline int32_t mixer_clamp(int32_t v, int32_t min, int32_t max) { return min > v ? min : (v > max ? max : v); }
static inline int mixer_min(int a, int b) { return (a < b) ? a : b; }
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
		tinymixer_free(buffer);
	}
}

static Source *find_source() {
	Source* best_source = 0;
	int32_t best_remaining = 0x7fffffff;

	for (int ii = 0; ii < c_nsources; ++ii) {
		Source* source = &g_mixer.sources[ii];
		if (!source->buffer)
			return source;

		if (0 == (source->flags & SourceFlags::Looping)) {
			const int32_t remaining = source->buffer->nsamples - source->sample_pos;
			if (remaining < best_remaining) {
				best_source = source;
				best_remaining = remaining;
			}
		}
	}

	if (NULL != best_source) {
		if (best_source->buffer) {
			decref((Buffer*)best_source->buffer);
			best_source->buffer = NULL;
		}
	}
	return best_source;
}

static int32_t apply_gain(int16_t sample, int qgain) {
	return (sample * qgain) / c_quantize;
}

static void render(Source* source, int32_t* buffer, const int qgain[2]) {
	const bool looping = 0 != (source->flags & SourceFlags::Looping);
	const int nchannels = source->buffer->nchannels;

	int remaining = c_nsamples;
	while (remaining) {
		int samples_read = mixer_min(source->buffer->nsamples - source->sample_pos, remaining);
		int samples_written = samples_read;
		const int16_t* samples = (const int16_t*)(source->buffer + 1) + (source->sample_pos * nchannels);

		// source has a non-1.0f frequency shift
		if (source->flags & SourceFlags::Frequency) {
			// samples_read/_written need to be modified to compensate for the frequency difference
			const int qfreq = (int)(source->frequency * c_quantize);
			samples_read = (samples_read * qfreq) / c_quantize;
			if (samples_read >= source->buffer->nsamples - source->sample_pos) {
				// adjust samples_written so we loop properly if needed
				samples_read = source->buffer->nsamples - source->sample_pos;
				samples_written = (samples_read * c_quantize) / qfreq;
				if (samples_read & c_quantize_mask) {
					++samples_written;
				}
			}

			for (int ii = 0; ii < samples_written; ++ii) {
				const int offset0 = ((ii+0) * qfreq) / c_quantize;
				int offset1 = ((ii+1) * qfreq) / c_quantize;
				// handle the case where the interpolation target is out of range (loop or clamp depending on source)
				if (offset1 >= source->buffer->nsamples - source->sample_pos) {
					offset1 = (looping ? 0 : source->buffer->nsamples - source->sample_pos - 1);
				}

				const int qinterp = (ii*qfreq) & c_quantize_mask;

				int32_t new_samples[2];
				if (nchannels == 1) {
					const int32_t sample = samples[offset0] + qinterp * (samples[offset1] - samples[offset0]) / c_quantize;
					new_samples[0] = new_samples[1] = sample;
				} else {
					new_samples[0] = samples[2*offset0 + 0] + qinterp * (samples[2*offset1 + 0] - samples[2*offset0 + 0]) / c_quantize;
					new_samples[1] = samples[2*offset0 + 1] + qinterp * (samples[2*offset1 + 1] - samples[2*offset0 + 1]) / c_quantize;
				}

				buffer[0] += apply_gain(new_samples[0], qgain[0]);
				buffer[1] += apply_gain(new_samples[1], qgain[1]);
				buffer += 2;
			}
		} else {
			for (int ii = 0; ii < samples_written; ++ii) {
				if (nchannels == 1) {
					buffer[0] += apply_gain(samples[ii], qgain[0]);
					buffer[1] += apply_gain(samples[ii], qgain[1]);
				} else {
					buffer[0] += apply_gain(samples[2*ii + 0], qgain[0]);
					buffer[1] += apply_gain(samples[2*ii + 1], qgain[1]);
				}

				buffer += 2;
			}
		}

		source->sample_pos += samples_read;
		remaining -= samples_written;

		if (remaining) {
			if (!looping)
				break;

			source->sample_pos = 0;
		}
	}
}

static void render_effects(int32_t* buffer) {

	int compressor_factor = g_mixer.compressor_factor;

	// get maximum absolute power level from the rendered buffer, and adjust the compressor factor
	int32_t max_power = 0;
	for (int ii = 0; ii < c_nsamples; ++ii) {
		const int32_t power = (int32_t)tinymixer_abs(buffer[ii]);
		if (power > max_power)
			max_power = power;
	}

	int target_compressor_factor = c_quantize;
	if (max_power > g_mixer.compressor_thresholds[1])
		target_compressor_factor = g_mixer.compressor_multipliers[1];
	else if (max_power > g_mixer.compressor_thresholds[0])
		target_compressor_factor = g_mixer.compressor_multipliers[0];

	int attack_release = c_quantize;
	if (target_compressor_factor < compressor_factor)
		attack_release = g_mixer.compressor_attack_per1ksamples;
	else if (target_compressor_factor > compressor_factor)
		attack_release = g_mixer.compressor_release_per1ksamples;

	// linearly interp compressor_factor toward the target compressor value
	const int32_t interp = (attack_release * c_nsamples) / (1000 * c_quantize);
	compressor_factor = compressor_factor + interp*(target_compressor_factor - compressor_factor);
	compressor_factor = mixer_clamp(compressor_factor, g_mixer.compressor_multipliers[1], c_quantize);

    // 2-pass compressor to limit dynamic range of audio clipping levels
	if (compressor_factor < c_quantize) {
		int32_t last_samples[2] = {
			g_mixer.compressor_last_samples[0],
			g_mixer.compressor_last_samples[1],
		};

		for (int ii = 0; ii < c_nsamples; ++ii) {
			const int32_t prev_sample = last_samples[ii&1];
			int32_t sample = buffer[ii];

			// uhhh... linear space? really??
			int32_t diff = sample - prev_sample;
			sample = prev_sample + (compressor_factor*diff)/ c_quantize;
			buffer[ii] = sample;

			last_samples[ii & 1] = sample;
		}


		g_mixer.compressor_last_samples[0] = last_samples[0];
		g_mixer.compressor_last_samples[1] = last_samples[1];
	}

	g_mixer.compressor_factor = compressor_factor;
}

static void mix(int32_t* buffer) {
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
	}

	tinymixer_memset(buffer, 0, sizeof(int32_t)*2*c_nsamples);

	// render playing sources
	for (int ii = 0; ii < nplaying; ++ii) {
		const int qgain[2] = {
			(int)(mixer_clamp(gain[ii][0], 0.0f, 1.0f) * c_quantize),
			(int)(mixer_clamp(gain[ii][1], 0.0f, 1.0f) * c_quantize),
		};

		Source* source = &g_mixer.sources[playing[ii]];
		render(source, buffer, qgain);
	}

	// allow application to apply a premixed track (such as music)
	if (g_mixer.callback)
		(g_mixer.callback)(buffer, c_nsamples, g_mixer.gain_callback);

	// render effects
	render_effects(buffer);

	// perform source-level post processing
	for (int ii = 0; ii < nplaying; ++ii) {
		Source* source = &g_mixer.sources[playing[ii]];

		// handle fadeout->stop
		if (source->flags & SourceFlags::FadeOut) {
			source->gain_base -= source->fadeout_per_sample * c_nsamples;
			if (source->gain_base <= 0.0f) {
				decref((Buffer*)source->buffer);
				source->buffer = 0;
				source->flags = 0;
			}
		} else if (source->sample_pos >= source->buffer->nsamples && 0 == (source->flags & SourceFlags::Looping)) {
			decref((Buffer*)source->buffer);
			source->buffer = 0;
			source->flags = 0;
		}
	}
}

void tinymixer_getsamples(int16_t* samples, int nsamples) {
	
	// was data leftover after the previous call to getsamples? Copy that out here
	while (nsamples && g_mixer.samples_remaining) {
		const int samples_to_mix = mixer_min(nsamples, g_mixer.samples_remaining);
		const int offset = 2*(c_nsamples - g_mixer.samples_remaining);

		// clip
		for (int ii = 0; ii < samples_to_mix*2; ++ii)
			samples[ii] = (int16_t)mixer_clamp(g_mixer.buffer[offset + ii], (int16_t)0x8000, 0x7fff);

		g_mixer.samples_remaining -= samples_to_mix;
		samples += (2*samples_to_mix);
		nsamples -= samples_to_mix;
	}

	// Copy out samples
	if (nsamples) {
		const int samples_to_mix = mixer_min(nsamples, c_nsamples);
		mix(g_mixer.buffer);
		g_mixer.samples_remaining = c_nsamples;

		// clip
		for (int ii = 0; ii < samples_to_mix*2; ++ii)
			samples[ii] = (int16_t)mixer_clamp(g_mixer.buffer[ii], (int16_t)0x8000, 0x7fff);


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
	source->sample_pos = 0;
	source->gain_base_index = (uint8_t)gain_index;

	addref((Buffer*)handle);
	return source;
}

static void play(Source* source) {
	source->flags |= SourceFlags::Playing;
}

void tinymixer_create_buffer_s16le(int channels, const int16_t* pcm_data, int pcm_data_size, const tinymixer_buffer** handle) {
	Buffer* buffer = (Buffer*)tinymixer_alloc(sizeof(Buffer) + pcm_data_size);
	buffer->refcnt = 1;
	buffer->nchannels = (uint8_t)channels;
	buffer->nsamples = pcm_data_size/sizeof(int16_t)/channels;

	// copy samples
	const int16_t* source = (const int16_t*)pcm_data;
	int16_t* dest = (int16_t*)(buffer + 1);
	const int nsamples = pcm_data_size/sizeof(uint16_t);
	for (int ii = 0; ii < nsamples; ++ii)
		*dest++ = *source++;

	*handle = (tinymixer_buffer*)buffer;
}

void tinymixer_create_buffer_float(int channels, const float* pcm_data, int pcm_data_size, const tinymixer_buffer** handle) {
	Buffer* buffer = (Buffer*)tinymixer_alloc(sizeof(Buffer) + pcm_data_size);
	buffer->refcnt = 1;
	buffer->nchannels = (uint8_t)channels;
	buffer->nsamples = pcm_data_size/sizeof(float)/channels;

	// copy samples
	const float* source = (const float*)pcm_data;
	int16_t* dest = (int16_t*)(buffer + 1);
	const int nsamples = pcm_data_size/sizeof(float);
	for (int ii = 0; ii < nsamples; ++ii)
		*dest++ = (int16_t)(*source++ * (float)0x8000);

	*handle = (tinymixer_buffer*)buffer;
}

int tinymixer_get_buffer_size(const tinymixer_buffer* handle) {
	const Buffer* buffer = (const Buffer*)handle;
	return sizeof(Buffer)+sizeof(int16_t)*buffer->nchannels*buffer->nsamples;
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
	source->flags &= ~SourceFlags::Looping;

	// set the current position to the end of the source and let mixing loop cleanup the source
	source->sample_pos = source->buffer->nsamples;
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
	g_mixer.compressor_factor = c_quantize;
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
	g_mixer.compressor_thresholds[0] = (int32_t)(0x8000 * mixer_clamp(thresholds[0], 0.0f, 1.0f));
	g_mixer.compressor_thresholds[1] = (int32_t)(0x8000 * mixer_clamp(thresholds[1], 0.0f, 1.0f));
	g_mixer.compressor_multipliers[0] = (int32_t)(c_quantize * mixer_clamp(multipliers[0], 0.0f, 1.0f));
	g_mixer.compressor_multipliers[1] = (int32_t)(c_quantize * mixer_clamp(multipliers[1], 0.0f, 1.0f));
	g_mixer.compressor_attack_per1ksamples = (int32_t)(1000.0f * c_quantize / (attack_seconds * g_mixer.sample_rate));
	g_mixer.compressor_release_per1ksamples = (int32_t)(1000.0f * c_quantize / (release_seconds * g_mixer.sample_rate));
}

void tinymixer_stop_all_sources() {
	for (int ii = 0; ii < c_nsources; ++ii) {
		Source* source = &g_mixer.sources[ii];
		if (source->buffer) {
			decref((Buffer*)source->buffer);
			source->buffer = 0;
			source->flags = 0;
		}
	}
}
