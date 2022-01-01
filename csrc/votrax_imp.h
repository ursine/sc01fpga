#pragma once

#include <stdint.h>
#include <stdbool.h>

// useful macros to deal with bit shuffling encryptions
#define BIT(x,n) (((x)>>(n))&1)

#define BITSWAP4(val,B3,B2,B1,B0) \
	((BIT(val,B3) << 3) | (BIT(val,B2) << 2) | (BIT(val,B1) << 1) | (BIT(val,B0) << 0))

#define BITSWAP7(val,B6,B5,B4,B3,B2,B1,B0) \
	((BIT(val,B6) << 6) | (BIT(val,B5) << 5) | (BIT(val,B4) << 4) | \
		(BIT(val,B3) << 3) | (BIT(val,B2) << 2) | (BIT(val,B1) << 1) | (BIT(val,B0) << 0))

#define CLEAR_LINE false
#define ASSERT_LINE true


typedef enum {
    T_COMMIT_PHONE,
    T_END_OF_PHONE
} clock_enum;

typedef struct mame_timer
{
    struct mame_timer *next;
    struct mame_timer *prev;
    void (*callback)(int);
    int callback_param;
    int tag;
    uint8_t enabled;
    uint8_t temporary;
    double period;
    double start;
    double expire;
} mame_timer;

typedef struct votraxsc01_vars {
    mame_timer* timer;

    uint32_t mainclock;                                // Current main clock
    double   sclock;                                   // Stream sample clock (40KHz, main/18)
    double   cclock;                                   // 20KHz capacitor switching clock (main/36)
    uint32_t sample_count;                             // Sample counter, to cadence chip updates

    // Inputs
    uint8_t  inflection;                                // 2-bit inflection value
    uint8_t  phone;                                     // 6-bit phone value

    // Outputs
    //!! devcb_write_line m_ar_cb;                      // Callback for ar
    bool     ar_state;                                  // Current ar state

    // "Unpacked" current rom values
    uint8_t  rom_duration;                              // Duration in 5KHz units (main/144) of one tick, 16 ticks per phone, 7 bits
    uint8_t  rom_vd, rom_cld;                           // Duration in ticks of the "voice" and "closure" delays, 4 bits
    uint8_t  rom_fa, rom_fc, rom_va;                    // Analog parameters, noise volume, noise freq cutoff and voice volume, 4 bits each
    uint8_t  rom_f1, rom_f2, rom_f2q, rom_f3;           // Analog parameters, formant frequencies and Q, 4 bits each
    bool     rom_closure;                               // Closure bit, true = silence at cld
    bool     rom_pause;                                 // Pause bit

    // Current interpolated values (8 bits each)
    uint8_t cur_fa, cur_fc, cur_va;
    uint8_t cur_f1, cur_f2, cur_f2q, cur_f3;

    // Current committed values
    uint8_t filt_fa, filt_fc, filt_va;                 // Analog parameters, noise volume, noise freq cutoff and voice volume, 4 bits each
    uint8_t filt_f1, filt_f2, filt_f2q, filt_f3;       // Analog parameters, formant frequencies/Q on 4 bits except f2 on 5 bits

    // Internal counters
    uint16_t phonetick;                                // 9-bits phone tick duration counter
    uint8_t  ticks;                                    // 5-bits tick counter
    uint8_t  pitch;                                    // 7-bits pitch counter
    uint8_t  closure;                                  // 5-bits glottal closure counter
    uint8_t  update_counter;                           // 6-bits counter for the 625Hz (main/1152) and 208Hz (main/3456) update timing generators

    // Internal state
    bool     cur_closure;                              // Current internal closure state
    uint16_t noise;                                    // 15-bit noise shift register
    bool     cur_noise;                                // Current noise output

    // Filter coefficients and level histories
    double voice_1[4];
    double voice_2[4];
    double voice_3[4];

    double noise_1[3];
    double noise_2[3];
    double noise_3[2];
    double noise_4[2];

    double vn_1[4];
    double vn_2[4];
    double vn_3[4];
    double vn_4[4];
    double vn_5[2];
    double vn_6[2];

    double f1_a[4],  f1_b[4];                   // F1 filtering
    double f2v_a[4], f2v_b[4];                  // F2 voice filtering
    double f2n_a[2], f2n_b[2];                  // F2 noise filtering
    double f3_a[4],  f3_b[4];                   // F3 filtering
    double f4_a[4],  f4_b[4];                   // F4 filtering
    double fx_a[1],  fx_b[2];                   // Final filtering
    double fn_a[3],  fn_b[3];                   // Noise shaping

    int stream;
} votraxsc01_vars;

extern int votrax_start();
void votraxsc01_w(uint8_t);

void Votrax_Update(int16_t*, size_t);

extern void timer_init();

// Filter updates
void filters_commit(int);

mame_timer* timer_alloc(void (*callback)(int));




#define TIME_IN_HZ(hz)        (1.0 / (double)(hz))
//#define TIME_IN_CYCLES(c,cpu) ((double)(c) * cycles_to_sec[cpu])
#define TIME_IN_SEC(s)        ((double)(s))
#define TIME_IN_MSEC(ms)      ((double)(ms) * (1.0 / 1000.0))
#define TIME_IN_USEC(us)      ((double)(us) * (1.0 / 1000000.0))
#define TIME_IN_NSEC(us)      ((double)(us) * (1.0 / 1000000000.0))

#define TIME_NOW              (0.0)
#define TIME_NEVER            (1.0e30)
