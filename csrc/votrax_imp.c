#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <time.h>
#include "votrax_imp.h"

// Filter Constants
static const unsigned int f1_caps[4]   = { 2546, 4973, 9861, 19724 };
static const unsigned int f2v1_caps[4] = { 1390, 2965, 5875, 11297 };
static const unsigned int f2v2_caps[5] = { 833, 1663, 3164, 6327, 12654 };
static const unsigned int f2n1_caps[4] = { 1390, 2965, 5875, 11297 };
static const unsigned int f2n2_caps[5] = { 833, 1663, 3164, 6327, 12654 };
static const unsigned int f3_caps[4]   = { 2226, 4485, 9056, 18111 };

static const double s_glottal_wave[] = {
    0,
    -4/7.0,
    7/7.0,
    6/7.0,
    5/7.0,
    4/7.0,
    3/7.0,
    2/7.0,
    1/7.0
};


// ROM definition for the Votrax phone ROM
static const uint8_t sc01a_bin[512] = {
        0xA4, 0x50, 0xA0, 0xF0, 0xE0, 0x00, 0x00, 0x03, 0xA4, 0x50, 0xA0, 0x00,
        0x23, 0x0A, 0x00, 0x3E, 0xA4, 0x58, 0xA0, 0x30, 0xF0, 0x00, 0x00, 0x3F,
        0xA3, 0x80, 0x69, 0xB0, 0xC1, 0x0C, 0x00, 0x3D, 0x26, 0xD3, 0x49, 0x90,
        0xA1, 0x09, 0x00, 0x3C, 0x27, 0x81, 0x68, 0x94, 0x21, 0x0A, 0x00, 0x3B,
        0x82, 0xC3, 0x48, 0x24, 0xA1, 0x08, 0x00, 0x3A, 0xA4, 0x00, 0x38, 0x18,
        0x68, 0x01, 0x00, 0x39, 0x20, 0x52, 0xE1, 0x88, 0x63, 0x0A, 0x00, 0x38,
        0x22, 0xC1, 0xE8, 0x90, 0x61, 0x04, 0x00, 0x37, 0xA2, 0x83, 0x60, 0x10,
        0x66, 0x03, 0x00, 0x36, 0xA2, 0xC1, 0xE8, 0x80, 0xA1, 0x09, 0x00, 0x35,
        0xA2, 0xC1, 0xE8, 0x34, 0x61, 0x0A, 0x00, 0x34, 0xA3, 0x81, 0x89, 0xB4,
        0x21, 0x0A, 0x00, 0x33, 0xA3, 0x81, 0x89, 0xE4, 0xA1, 0x07, 0x00, 0x32,
        0xA3, 0x81, 0x89, 0x54, 0x63, 0x01, 0x00, 0x31, 0xA3, 0x80, 0x69, 0x60,
        0x61, 0x04, 0x00, 0x30, 0xA7, 0x80, 0xE8, 0x74, 0xA0, 0x07, 0x00, 0x2F,
        0xA7, 0x80, 0xE8, 0x74, 0x20, 0x0A, 0x00, 0x2E, 0x22, 0xC1, 0x60, 0x14,
        0x66, 0x0A, 0x00, 0x2D, 0x26, 0xD3, 0x49, 0x70, 0x20, 0x0A, 0x00, 0x2C,
        0x82, 0x43, 0x08, 0x54, 0x63, 0x04, 0x00, 0x2B, 0xE0, 0x32, 0x11, 0xE8,
        0x72, 0x01, 0x00, 0x2A, 0x26, 0x53, 0x01, 0x64, 0xA1, 0x07, 0x00, 0x29,
        0x22, 0xC1, 0xE8, 0x80, 0x21, 0x0A, 0x00, 0x28, 0xA6, 0x91, 0x61, 0x80,
        0x21, 0x0A, 0x00, 0x27, 0xA2, 0xC1, 0xE8, 0x84, 0x21, 0x0A, 0x00, 0x26,
        0xA8, 0x24, 0x13, 0x63, 0xB2, 0x07, 0x00, 0x25, 0xA3, 0x40, 0xE9, 0x84,
        0xC1, 0x0C, 0x00, 0x24, 0xA3, 0x81, 0x89, 0x54, 0xE3, 0x00, 0x00, 0x23,
        0x26, 0x12, 0xA0, 0x64, 0x61, 0x0A, 0x00, 0x22, 0x26, 0xD3, 0x69, 0x70,
        0x61, 0x05, 0x00, 0x21, 0xA6, 0xC1, 0xC9, 0x84, 0x21, 0x0A, 0x00, 0x20,
        0xE0, 0x32, 0x91, 0x48, 0x68, 0x04, 0x00, 0x1F, 0x26, 0x91, 0xE8, 0x00,
        0x7C, 0x0B, 0x00, 0x1E, 0xA8, 0x2C, 0x83, 0x65, 0xA2, 0x07, 0x00, 0x1D,
        0x26, 0xC1, 0x41, 0xE0, 0x73, 0x01, 0x00, 0x1C, 0xAC, 0x04, 0x22, 0xFD,
        0x62, 0x01, 0x00, 0x1B, 0x2C, 0x34, 0x7B, 0xDB, 0xE8, 0x00, 0x00, 0x1A,
        0x2C, 0x64, 0x23, 0x11, 0x72, 0x0A, 0x00, 0x19, 0xA2, 0xD0, 0x09, 0xF4,
        0xA1, 0x07, 0x00, 0x18, 0x23, 0x81, 0x49, 0x20, 0x21, 0x0A, 0x00, 0x17,
        0x23, 0x81, 0x49, 0x30, 0xA1, 0x07, 0x00, 0x16, 0xA3, 0x40, 0xE9, 0x84,
        0xA1, 0x08, 0x00, 0x15, 0x36, 0x4B, 0x08, 0xD4, 0xA0, 0x09, 0x00, 0x14,
        0xA3, 0x80, 0x69, 0x70, 0xA0, 0x08, 0x00, 0x13, 0x60, 0x58, 0xD1, 0x9C,
        0x63, 0x01, 0x00, 0x12, 0x6C, 0x54, 0x8B, 0xFB, 0xA2, 0x09, 0x00, 0x11,
        0x6C, 0x54, 0x8B, 0xFB, 0x63, 0x01, 0x00, 0x10, 0x28, 0x64, 0xD3, 0xF7,
        0x63, 0x01, 0x00, 0x0F, 0x22, 0x91, 0xE1, 0x90, 0x73, 0x01, 0x00, 0x0E,
        0x36, 0x19, 0x24, 0xE6, 0x61, 0x0A, 0x00, 0x0D, 0x32, 0x88, 0xA5, 0x66,
        0xA3, 0x07, 0x00, 0x0C, 0xA6, 0x91, 0x61, 0x90, 0xA1, 0x09, 0x00, 0x0B,
        0xA6, 0x91, 0x61, 0x90, 0x61, 0x0A, 0x00, 0x0A, 0xA6, 0x91, 0x61, 0x80,
        0x61, 0x0B, 0x00, 0x09, 0xA3, 0x40, 0xE9, 0xC4, 0x61, 0x01, 0x00, 0x08,
        0x6C, 0x54, 0xCB, 0xF3, 0x63, 0x04, 0x00, 0x07, 0xA6, 0xC1, 0xC9, 0x34,
        0xA1, 0x07, 0x00, 0x06, 0xA6, 0xC1, 0xC9, 0x64, 0x61, 0x01, 0x00, 0x05,
        0xE8, 0x16, 0x03, 0x61, 0xFB, 0x00, 0x00, 0x04, 0x27, 0x81, 0x68, 0xC4,
        0xA1, 0x09, 0x00, 0x02, 0x27, 0x81, 0x68, 0xD4, 0x61, 0x01, 0x00, 0x01,
        0x27, 0x81, 0x68, 0x74, 0x61, 0x03, 0x00, 0x00
};

static const char *PhonemeNames[65] = {
    "EH3",  "EH2",  "EH1",  "PA0",  "DT",   "A1",   "A2",   "ZH",
    "AH2",  "I3",   "I2",   "I1",   "M",    "N",    "B",    "V",
    "CH",   "SH",   "Z",    "AW1",  "NG",   "AH1",  "OO1",  "OO",
    "L",    "K",    "J",    "H",    "G",    "F",    "D",    "S",
    "A",    "AY",   "Y1",   "UH3",  "AH",   "P",    "O",    "I",
    "U",    "Y",    "T",    "R",    "E",    "W",    "AE",   "AE1",
    "AW2",  "UH2",  "UH1",  "UH",   "O2",   "O1",   "IU",   "U1",
    "THV",  "TH",   "ER",   "EH",   "E1",   "AW",   "PA1",  "STOP",
    0
};

// Compute a total capacitor value based on which bits are currently active
static unsigned int bits_to_caps(unsigned int value, const unsigned int* const caps_values, const size_t N) {
    size_t i;
    unsigned int total = 0;
    for(i = 0; i < N; ++i) {
        if(value & 1)
            total += caps_values[i];
        value >>= 1;
    }
    return total;
}


static struct votraxsc01_vars votraxsc01_locals;

/*
  Playing with analog filters, or where all the magic filter formulas are coming from.

  First you start with an analog circuit, for instance this one:

  |                     +--[R2]--+
  |                     |        |
  |                     +--|C2|--+<V1     +--|C3|--+
  |                     |        |        |        |
  |  Vi   +--[R1]--+    |  |\    |        |  |\    |
  |  -----+        +----+--+-\   |        +--+-\   |
  |       +--|C1|--+       |  >--+--[Rx]--+  |  >--+----- Vo
  |                |     0-++/             0-++/   |
  |                |       |/    +--[R0]--+  |/    |
  |                |             |        |        |
  |                |             |    /|  |        |
  |                |             |   /-+--+--[R0]--+
  |                +--[R4]-------+--<  |
  |                            V2^   \++-0
  |                                   \|

  It happens to be what most of the filters in the sc01a look like.

  You need to determine the transfer function H(s) of the circuit, which is
  defined as the ratio Vo/Vi.  To do that, you use some properties:

  - The intensity through an element is equal to the voltage
    difference through the element divided by the impedence

  - The impedence of a resistance is equal to its resistance

  - The impedence of a capacitor is 1/(s*C) where C is its capacitance

  - The impedence of elements in series is the sum of the impedences

  - The impedence of elements in parallel is the inverse of the sum of
    the inverses

  - The sum of all intensities flowing into a node is 0 (there's no
    charge accumulation in a wire)

  - An operational amplifier in looped mode is an interesting beast:
    the intensity at its two inputs is always 0, and the voltage is
    forced identical between the inputs.  In our case, since the '+'
    inputs are all tied to ground, that means that the '-' inputs are at
    voltage 0, intensity 0.

  From here we can build some equations.  Noting:
  X1 = 1/(1/R1 + s*C1)
  X2 = 1/(1/R2 + s*C2)
  X3 = 1/(s*C3)

  Then computing the intensity flow at each '-' input we have:
  Vi/X1 + V2/R4 + V1/X2 = 0
  V2/R0 + Vo/R0 = 0
  V1/Rx + Vo/X3 = 0

  Wrangling the equations, one eventually gets:
  |                            1 + s * C1*R1
  | Vo/Vi = H(s) = (R4/R1) * -------------------------------------------
  |                            1 + s * C3*Rx*R4/R2 + s^2 * C2*C3*Rx*R4

  To check the mathematics between the 's' stuff, check "Laplace
  transform".  In short, it's a nice way of manipulating derivatives
  and integrals without having to manipulate derivatives and
  integrals.

  With that transfer function, we first can compute what happens to
  every frequency in the input signal.  You just compute H(2i*pi*f)
  where f is the frequency, which will give you a complex number
  representing the amplitude and phase effect.  To get the usual dB
  curves, compute 20*log10(abs(v))).

  Now, once you have an analog transfer function, you can build a
  digital filter from it using what is called the bilinear transform.

  In our case, we have an analog filter with the transfer function:
  |                 1 + k[0]*s
  |        H(s) = -------------------------
  |                 1 + k[1]*s + k[2]*s^2

  We can always reintroduce the global multipler later, and it's 1 in
  most of our cases anyway.

  The we pose:
  |                    z-1
  |        s(z) = zc * ---
  |                    z+1

  where zc = 2*pi*fr/tan(pi*fr/fs)
  with fs = sampling frequency
  and fr = most interesting frequency

  Then we rewrite H in function of negative integer powers of z.

  Noting m0 = zc*k[0], m1 = zc*k[1], m2=zc*zc*k[2],

  a little equation wrangling then gives:

  |                 (1+m0)    + (3+m0)   *z^-1 + (3-m0)   *z^-2 +    (1-m0)*z^-3
  |        H(z) = ----------------------------------------------------------------
  |                 (1+m1+m2) + (3+m1-m2)*z^-1 + (3-m1-m2)*z^-2 + (1-m1+m2)*z^-3

  That beast in the digital transfer function, of which you can
  extract response curves by posing z = exp(2*i*pi*f/fs).

  Note that the bilinear transform is an approximation, and H(z(f)) =
  H(s(f)) only at frequency fr.  And the shape of the filter will be
  better respected around fr.  If you look at the curves of the
  filters we're interested in, the frequency:
  fr = sqrt(abs(k[0]*k[1]-k[2]))/(2*pi*k[2])

  which is a (good) approximation of the filter peak position is a
  good choice.

  Note that terminology wise, the "standard" bilinear transform is
  with fr = fs/2, and using a different fr is called "pre-warping".

  So now we have a digital transfer function of the generic form:

  |                 a[0] + a[1]*z^-1 + a[2]*z^-2 + a[3]*z^-3
  |        H(z) = --------------------------------------------
  |                 b[0] + b[1]*z^-1 + b[2]*z^-2 + b[3]*z^-3

  The magic then is that the powers of z represent time in samples.
  Noting x the input stream and y the output stream, you have:
  H(z) = y(z)/x(z)

  or in other words:
  y*b[0]*z^0 + y*b[1]*z^-1 + y*b[2]*z^-2 + y*b[3]*z^-3 = x*a[0]*z^0 + x*a[1]*z^-1 + x*a[2]*z^-2 + x*a[3]*z^-3

  i.e.

  y*z^0 = (x*a[0]*z^0 + x*a[1]*z^-1 + x*a[2]*z^-2 + x*a[3]*z^-3 - y*b[1]*z^-1 - y*b[2]*z^-2 - y*b[3]*z^-3) / b[0]

  and powers of z being time in samples,

  y[0] = (x[0]*a[0] + x[-1]*a[1] + x[-2]*a[2] + x[-3]*a[3] - y[-1]*b[1] - y[-2]*b[2] - y[-3]*b[3]) / b[0]

  So you have a filter you can apply.  Note that this is why you want
  negative powers of z.  Positive powers would mean looking into the
  future (which is possible in some cases, in particular with x, and
  has some very interesting properties, but is not very useful in
  analog circuit simulation).

  Note that if you have multiple inputs, all this stuff is linear.
  Or, in other words, you just have to split it in multiple circuits
  with only one input connected each time and sum the results.  It
  will be correct.

  Also, since we're in practice in a dynamic system, for an amplifying
  filter (i.e. where things like r4/r1 is not 1), it's better to
  proceed in two steps:

  - amplify the input by the current value of the coefficient, and
    historize it
  - apply the now non-amplifying filter to the historized amplified
    input

  That way reduces the probability of the output bouncing all over the
  place.

  Except, we're not done yet.  Doing resistors precisely in an IC is
  very hard and/or expensive (you may have heard of "laser cut
  resistors" in DACs of the time).  Doing capacitors is easier, and
  their value is proportional to their surface.  So there are no
  resistors on the sc01 die (which is a lie, there are three, but not
  in the filter path.  They are used to scale the voltage in the pitch
  wave and to generate +5V from the +9V), but a magic thing called a
  switched capacitor.  Lookup patent 4,433,210 for details.  Using
  high frequency switching a capacitor can be turned into a resistor
  of value 1/(C*f) where f is the switching frequency (20Khz,
  main/36).  And the circuit is such that the absolute value of the
  capacitors is irrelevant, only their ratio is useful, which factors
  out the intrinsic capacity-per-surface-area of the IC which may be
  hard to keep stable from one die to another.  As a result all the
  capacitor values we use are actually surfaces in square micrometers.

  For the curious, it looks like the actual capacitance was around 25
  femtofarad per square micrometer.

*/

static void build_standard_filter(double* const a,
                                  double* const b,
                                  const double c1t, // Unswitched cap, input, top
                                  const double c1b, // Switched cap, input, bottom
                                  const double c2t, // Unswitched cap, over first amp-op, top
                                  const double c2b, // Switched cap, over first amp-op, bottom
                                  const double c3,  // Cap between the two op-amps
                                  const double c4)  // Cap over second op-amp
{
    // First compute the three coefficients of H(s).  One can note
    // that there is as many capacitor values on both sides of the
    // division, which confirms that the capacity-per-surface-area
    // is not needed.
    const double k0 = c1t / (votraxsc01_locals.cclock * c1b);
    const double k1 = c4 * c2t / (votraxsc01_locals.cclock * c1b * c3);
    const double k2 = c4 * c2b / (votraxsc01_locals.cclock * votraxsc01_locals.cclock * c1b * c3);

    // Estimate the filter cutoff frequency
    const double fpeak = sqrt(fabs(k0*k1 - k2))/((2.*M_PI)*k2);

    // figure radians from frequency
    double w0 = fpeak * (2.0 * M_PI);

    // keep it in -PI/T .. PI/T
    const double wdMax = M_PI * votraxsc01_locals.sclock;
    if (w0 > wdMax) w0 -= 2.0 * wdMax;

    // Turn that into a warp multiplier
    const double zc = w0/tan(w0 / (2.0*votraxsc01_locals.sclock));

    // Finally compute the result of the z-transform
    const double m0 = zc*k0;
    const double m1 = zc*k1;
    const double m2 = zc*zc*k2;

    a[0] = 1.+m0;
    a[1] = 3.+m0;
    a[2] = 3.-m0;
    a[3] = 1.-m0;
    b[0] = 1.+m1+m2;
    b[1] = 3.+m1-m2;
    b[2] = 3.-m1-m2;
    b[3] = 1.-m1+m2;
}


/*
  Second filter type used once at the end, much simpler:

  |           +--[R1]--+
  |           |        |
  |           +--|C1|--+
  |           |        |
  |  Vi       |  |\    |
  |  ---[R0]--+--+-\   |
  |              |  >--+------ Vo
  |            0-++/
  |              |/


  Vi/R0 = Vo / (1/(1/R1 + s.C1)) = Vo (1/R1 + s.C1)
  H(s) = Vo/Vi = (R1/R0) * (1 / (1 + s.R1.C1))
*/

static void build_lowpass_filter(double * const a, double * const b,
                                 const double c1t, // Unswitched cap, over amp-op, top
                                 const double c1b) // Switched cap, over amp-op, bottom
{
    // The cap value puts the cutoff at around 150Hz, put that's no good.
    // Recordings show we want it around 4K, so fuzz it.

    // Compute the only coefficient we care about
    const double k = c1b / (votraxsc01_locals.cclock * c1t) * (150.0/4000.0);

    // Compute the filter cutoff frequency
    const double fpeak = 1./((2.*M_PI)*k);

    // figure radians from frequency
    double w0 = fpeak * (2.0 * M_PI);

    // keep it in -PI/T .. PI/T
    const double wdMax = M_PI * votraxsc01_locals.sclock;
    if (w0 > wdMax)
        w0 -= 2.0 * wdMax;

    // Turn that into a warp multiplier
    const double zc = w0/tan(w0 / (2.0*votraxsc01_locals.sclock));

    // Finally, compute the result of the z-transform
    const double m = zc*k;

    a[0] = 1.;
    b[0] = 1.+m;
    b[1] = 1.-m;
}

/*
  Used to shape the white noise

         +-------------------------------------------------------------------+
         |                                                                   |
         +--|C1|--+---------|C3|----------+--|C4|--+                         |
         |        |      +        +       |        |                         |
   Vi    |  |\    |     (1)      (1)      |        |       +        +        |
   -|R0|-+--+-\   |      |        |       |  |\    |      (1)      (1)       |
            |  >--+--(2)-+--|C2|--+---(2)-+--+-\   |       |        |        |
          0-++/          |                   |  >--+--(2)--+--|C5|--+---(2)--+
            |/          Vo                 0-++/
                                             |/
   Equivalent:

         +------------------|R5|-------------------+
         |                                         |
         +--|C1|--+---------|C3|----------+--|C4|--+
         |        |                       |        |
   Vi    |  |\    |                       |        |
   -|R0|-+--+-\   |                       |  |\    |
            |  >--+---------|R2|----------+--+-\   |
          0-++/   |                          |  >--+
            |/   Vo                        0-++/
                                             |/

  We assume r0 = r2
*/

static void build_noise_shaper_filter(double* const a,
                                      double* const b,
                                      const double c1,  // Cap over first amp-op
                                      const double c2t, // Unswitched cap between amp-ops, input, top
                                      const double c2b, // Switched cap between amp-ops, input, bottom
                                      const double c3,  // Cap over second amp-op
                                      const double c4)  // Switched cap after second amp-op
{
    // Coefficients of H(s) = k1*s / (1 + k2*s + k3*s^2)
    const double k0 = c2t*c3*c2b/c4;
    const double k1 = c2t*(votraxsc01_locals.cclock * c2b);
    const double k2 = c1*c2t*c3/(votraxsc01_locals.cclock * c4);

    // Estimate the filter cutoff frequency
    const double fpeak = sqrt(1./k2)/(2.*M_PI);

    // figure radians from frequency
    double w0 = fpeak * (2.0 * M_PI);

    // keep it in -PI/T .. PI/T
    const double wdMax = M_PI * votraxsc01_locals.sclock;
    if (w0 > wdMax) w0 -= 2.0 * wdMax;

    // Turn that into a warp multiplier
    const double zc = w0/tan(w0 / (2.0*votraxsc01_locals.sclock));

    // Finally compute the result of the z-transform
    const double m0 = zc*k0;
    const double m1 = zc*k1;
    const double m2 = zc*zc*k2;

    a[0] = m0;
    a[1] = 0;
    a[2] = -m0;
    b[0] = 1.+m1+m2;
    b[1] = 2.-2.*m2;
    b[2] = 1.-m1+m2;
}

/*
  Noise injection in f2

  |                     +--[R2]--+        +--[R1]-------- Vi
  |                     |        |        |
  |                     +--|C2|--+<V1     +--|C3|--+
  |                     |        |        |        |
  |                     |  |\    |        |  |\    |
  |                +----+--+-\   |        +--+-\   |
  |                |       |  >--+--[Rx]--+  |  >--+----- Vo
  |                |     0-++/             0-++/   |
  |                |       |/    +--[R0]--+  |/    |
  |                |             |        |        |
  |                |             |    /|  |        |
  |                |             |   /-+--+--[R0]--+
  |                +--[R4]-------+--<  |
  |                            V2^   \++-0
  |                                   \|

  We drop r0/r1 out of the equation (it factorizes), and we rescale so
  that H(infinity)=1.
*/

static void build_injection_filter(double* const a,
                                   double* const b,
                                   const double c1b, // Switched cap, input, bottom
                                   const double c2t, // Unswitched cap, over first amp-op, top
                                   const double c2b, // Switched cap, over first amp-op, bottom
                                   const double c3,  // Cap between the two op-amps
                                   const double c4)  // Cap over second op-amp
{
    // First compute the three coefficients of H(s) = (k0 + k2*s)/(k1 - k2*s)
    const double k0 = votraxsc01_locals.cclock * c2t;
    const double k1 = votraxsc01_locals.cclock * (c1b * c3 / c2t - c2t);
    const double k2 = c2b;

    // Don't pre-warp
    const double zc = 2.*votraxsc01_locals.sclock;

    // Finally compute the result of the z-transform
    const double m = zc*k2;

    a[0] = k0 + m;
    a[1] = k0 - m;
    b[0] = k1 - m;
    b[1] = k1 + m;

    // That ends up in a numerically unstable filter.  Neutralize it for now.
    a[0] = 0;
    a[1] = 0;
    b[0] = 1.;
    b[1] = 0;
}

static void phone_commit()
{
    // Only these two counters are reset on phone change, the rest is
    // free-running.
    votraxsc01_locals.phonetick = 0;
    votraxsc01_locals.ticks = 0;

    // In the real chip, the rom is re-read all the time.  Since it's
    // internal and immutable, no point in not caching it though.
    for (int i = 0; i<64; i++) {
        const uint64_t val = ((uint64_t*)sc01a_bin)[i];
        if (votraxsc01_locals.phone == ((val >> 56) & 0x3f))
        {
            votraxsc01_locals.rom_f1  = BITSWAP4(val, 0, 7, 14, 21);
            votraxsc01_locals.rom_va  = BITSWAP4(val, 1, 8, 15, 22);
            votraxsc01_locals.rom_f2  = BITSWAP4(val, 2, 9, 16, 23);
            votraxsc01_locals.rom_fc  = BITSWAP4(val, 3, 10, 17, 24);
            votraxsc01_locals.rom_f2q = BITSWAP4(val, 4, 11, 18, 25);
            votraxsc01_locals.rom_f3  = BITSWAP4(val, 5, 12, 19, 26);
            votraxsc01_locals.rom_fa  = BITSWAP4(val, 6, 13, 20, 27);

            // These two values have their bit orders inverted
            // compared to everything else due to a bug in the
            // prototype (miswiring of the comparator with the ticks
            // count) they compensated in the rom.

            votraxsc01_locals.rom_cld = BITSWAP4(val, 34, 32, 30, 28);
            votraxsc01_locals.rom_vd  = BITSWAP4(val, 35, 33, 31, 29);

            votraxsc01_locals.rom_closure  = BIT(val, 36);
            votraxsc01_locals.rom_duration = BITSWAP7(~val, 37, 38, 39, 40, 41, 42, 43);

            // Hard-wired on the die, not an actual part of the rom.
            votraxsc01_locals.rom_pause = (votraxsc01_locals.phone == 0x03) || (votraxsc01_locals.phone == 0x3e);

            printf("commit fa=%x va=%x fc=%x f1=%x f2=%x f2q=%x f3=%x dur=%02x cld=%x vd=%d cl=%d pause=%d\n", votraxsc01_locals.rom_fa, votraxsc01_locals.rom_va, votraxsc01_locals.rom_fc, votraxsc01_locals.rom_f1, votraxsc01_locals.rom_f2, votraxsc01_locals.rom_f2q, votraxsc01_locals.rom_f3, votraxsc01_locals.rom_duration, votraxsc01_locals.rom_cld, votraxsc01_locals.rom_vd, votraxsc01_locals.rom_closure, votraxsc01_locals.rom_pause);

            // That does not happen in the sc01(a) rom, but let's
            // cover our behind.
            if (votraxsc01_locals.rom_cld == 0)
                votraxsc01_locals.cur_closure = votraxsc01_locals.rom_closure;

            return;
        }
    }
}

void filters_commit(int force)
{
    votraxsc01_locals.filt_fa = votraxsc01_locals.cur_fa >> 4;
    votraxsc01_locals.filt_fc = votraxsc01_locals.cur_fc >> 4;
    votraxsc01_locals.filt_va = votraxsc01_locals.cur_va >> 4;

    if (force || votraxsc01_locals.filt_f1 != votraxsc01_locals.cur_f1 >> 4) {
        votraxsc01_locals.filt_f1 = votraxsc01_locals.cur_f1 >> 4;

        build_standard_filter(votraxsc01_locals.f1_a, votraxsc01_locals.f1_b,
                              11247,
                              11797,
                              949,
                              52067,
                              2280 + bits_to_caps(votraxsc01_locals.filt_f1, f1_caps, 4),
                              166272);
    }

    if (force || votraxsc01_locals.filt_f2 != votraxsc01_locals.cur_f2 >> 3 || votraxsc01_locals.filt_f2q != votraxsc01_locals.cur_f2q >> 4) {
        votraxsc01_locals.filt_f2 = votraxsc01_locals.cur_f2 >> 3;
        votraxsc01_locals.filt_f2q = votraxsc01_locals.cur_f2q >> 4;

        build_standard_filter(votraxsc01_locals.f2v_a, votraxsc01_locals.f2v_b,
                              24840,
                              29154,
                              829 + bits_to_caps(votraxsc01_locals.filt_f2q, f2v1_caps, 4),
                              38180,
                              2352 + bits_to_caps(votraxsc01_locals.filt_f2, f2v2_caps, 5),
                              34270);

        build_injection_filter(votraxsc01_locals.f2n_a, votraxsc01_locals.f2n_b,
                               29154,
                               829 + bits_to_caps(votraxsc01_locals.filt_f2q, f2n1_caps, 4),
                               38180,
                               2352 + bits_to_caps(votraxsc01_locals.filt_f2, f2n2_caps, 5),
                               34270);
    }

    if (force || votraxsc01_locals.filt_f3 != votraxsc01_locals.cur_f3 >> 4) {
        votraxsc01_locals.filt_f3 = votraxsc01_locals.cur_f3 >> 4;
        build_standard_filter(votraxsc01_locals.f3_a, votraxsc01_locals.f3_b,
                              0,
                              17594,
                              868,
                              18828,
                              8480 + bits_to_caps(votraxsc01_locals.filt_f3, f3_caps, 4),
                              50019);
    }

    if (force) {
        build_standard_filter(votraxsc01_locals.f4_a, votraxsc01_locals.f4_b,
                              0,
                              28810,
                              1165,
                              21457,
                              8558,
                              7289);

        build_lowpass_filter(votraxsc01_locals.fx_a, votraxsc01_locals.fx_b,
                             1122,
                             23131);

        build_noise_shaper_filter(votraxsc01_locals.fn_a, votraxsc01_locals.fn_b,
                                  15500,
                                  14854,
                                  8450,
                                  9523,
                                  14083);
    }

    //if(votraxsc01_locals.filt_fa || votraxsc01_locals.filt_va || votraxsc01_locals.filt_fc || votraxsc01_locals.filt_f1 || votraxsc01_locals.filt_f2 || votraxsc01_locals.filt_f2q || votraxsc01_locals.filt_f3) {
        printf("filter fa=%x va=%x fc=%x f1=%x f2=%02x f2q=%x f3=%x\n",
                votraxsc01_locals.filt_fa, votraxsc01_locals.filt_va, votraxsc01_locals.filt_fc,
                votraxsc01_locals.filt_f1, votraxsc01_locals.filt_f2, votraxsc01_locals.filt_f2q,
                votraxsc01_locals.filt_f3);
    //}
}

int votrax_start()
{
    memset(&votraxsc01_locals, 0x00, sizeof(votraxsc01_locals));

    // Defines the sound interface
    //votraxsc01_locals.intf = msound->sound_interface;

    // initialize internal state
    // TODO: Figure out what a good starting clock is
    votraxsc01_locals.mainclock = 9000;
    votraxsc01_locals.sclock = votraxsc01_locals.mainclock / 18.0;
    votraxsc01_locals.cclock = votraxsc01_locals.mainclock / 36.0;

    // This is used to initialize the device
    //votraxsc01_locals.stream = stream_init_float("Votrax - SC01", votraxsc01_locals.intf->mixing_level[0], votraxsc01_locals.sclock, 0, Votrax_Update, 1);
    //votraxsc01_locals.timer = timer_alloc(VOTRAXSC01_sh_start_timeout);

    // reset outputs
    //!! m_ar_cb.resolve_safe();
    votraxsc01_locals.ar_state = ASSERT_LINE;

    //!! was the separate reset code from here on:

    // Technically, there's no reset in this chip, and initial state
    // is random.  Still, it's a good idea to start it with something
    // sane.

    votraxsc01_locals.phone = 0x3f;
    votraxsc01_locals.inflection = 0;
    votraxsc01_locals.ar_state = ASSERT_LINE;
    //!! m_ar_cb(votraxsc01_locals.ar_state);

    votraxsc01_locals.sample_count = 0;

    // Initialize the m_rom* values
    phone_commit();

    // Clear the interpolation sram
    votraxsc01_locals.cur_fa = votraxsc01_locals.cur_fc = votraxsc01_locals.cur_va = 0;
    votraxsc01_locals.cur_f1 = votraxsc01_locals.cur_f2 = votraxsc01_locals.cur_f2q = votraxsc01_locals.cur_f3 = 0;

    // Initialize the m_filt* values and the filter coefficients
    filters_commit(1);

    // Clear the rest of the internal digital state
    votraxsc01_locals.pitch = 0;
    votraxsc01_locals.closure = 0;
    votraxsc01_locals.update_counter = 0;
    votraxsc01_locals.cur_closure = 1;
    votraxsc01_locals.noise = 0;
    votraxsc01_locals.cur_noise = 0;

    return 0;
}

void votrax_stop()
{
    //if (votraxsc01_locals.timer )
    //    timer_remove(votraxsc01_locals.timer);
    //votraxsc01_locals.timer = 0;
}


void stream_update(int channel, int min_interval) {
    printf("Stream Update\n");
} // todo: output sound

void votraxsc01_w(const uint8_t data)
{
    uint8_t prev;

    // only 2 bits matter
    int inflection = (data >> 6) & 0x03; //!! MAME astrocde also uses: (data & 0x80) ? 0 : 2;
    //if (votraxsc01_locals.inflection == inflection) //!! original code, due to separate w for inflection in MAME -> does not matter, astrocde has inflection, then -directly- phone write after it
    //	return;

    stream_update(votraxsc01_locals.stream, 0);
    votraxsc01_locals.inflection = inflection;

    //!! in the original code this was the separate phone write from here on:

    // flush out anything currently processing
    //stream_update(votraxsc01_locals.stream,0);

    prev = votraxsc01_locals.phone;

    // only 6 bits matter
    votraxsc01_locals.phone = data & 0x3f;

    if (votraxsc01_locals.phone != prev || votraxsc01_locals.phone != 0x3f) {
        printf("phone %02x.%d %s\n", votraxsc01_locals.phone, votraxsc01_locals.inflection,
               PhonemeNames[votraxsc01_locals.phone]);
    }
    votraxsc01_locals.ar_state = CLEAR_LINE;

    // Schedule a commit/ar reset at roughly 0.1ms in the future (one
    // phi1 transition followed by the rom extra state in practice),
    // but only if there isn't already one on the fly.  It will
    // override an end-of-phone timeout if there's one pending, but
    // that's not a problem since stb does that anyway.
//    if (timer_expire(votraxsc01_locals.timer) == TIME_NEVER || timer_param(votraxsc01_locals.timer) != T_COMMIT_PHONE) {
//        timer_adjust(votraxsc01_locals.timer,
//                    72./(double)votraxsc01_locals.mainclock  /*attotime::from_ticks(72, votraxsc01_locals.mainclock)*/,
//                    T_COMMIT_PHONE, TIME_NEVER); //!! correct?
//    }
}


static void VOTRAXSC01_sh_start_timeout(const clock_enum which)
{
    stream_update(votraxsc01_locals.stream, 0);

    switch (which) {
        case T_COMMIT_PHONE:
            // strobe -> commit transition,
            phone_commit();
            //timer_adjust(votraxsc01_locals.timer, (double)(16 * (votraxsc01_locals.rom_duration * 4u + 1) * 4 * 9 + 2)/(double)votraxsc01_locals.mainclock/*attotime::from_ticks(16 * (votraxsc01_locals.rom_duration * 4 + 1) * 4 * 9 + 2, votraxsc01_locals.mainclock)*/, T_END_OF_PHONE, TIME_NEVER);
            break;

        case T_END_OF_PHONE:
            // end of phone
            votraxsc01_locals.ar_state = ASSERT_LINE;
            break;

        default:
            break;
    }

    //!! m_ar_cb(votraxsc01_locals.ar_state);
    //if (votraxsc01_locals.intf->BusyCallback[0])
    //    (*votraxsc01_locals.intf->BusyCallback[0])(!votraxsc01_locals.ar_state); //!! inverted behavior from MAME
}


static void interpolate(uint8_t* const reg, const uint8_t target)
{
    // One step of interpolation, adds one eight of the distance
    // between the current value and the target.
    *reg = *reg - (*reg >> 3) + (target << 1);
}


static void chip_update()
{
    // Phone tick counter update.  Stopped when ticks reach 16.
    // Technically the counter keeps updating, but the comparator is
    // disabled.
    if (votraxsc01_locals.ticks != 0x10) {
        votraxsc01_locals.phonetick++;
        // Comparator is with duration << 2, but there's a one-tick
        // delay in the path.
        if (votraxsc01_locals.phonetick == ((votraxsc01_locals.rom_duration << 2) | 1)) {
            votraxsc01_locals.phonetick = 0;
            votraxsc01_locals.ticks++;
            if (votraxsc01_locals.ticks == votraxsc01_locals.rom_cld)
                votraxsc01_locals.cur_closure = votraxsc01_locals.rom_closure;
        }
    }

    // The two update timing counters.  One divides by 16, the other
    // by 48, and they're phased so that the 208Hz counter ticks
    // exactly between two 625Hz ticks.
    votraxsc01_locals.update_counter++;
    if (votraxsc01_locals.update_counter == 0x30) {
        votraxsc01_locals.update_counter = 0;
    }

    int tick_625 = !(votraxsc01_locals.update_counter & 0xf);
    int tick_208 = (votraxsc01_locals.update_counter == 0x28);

    // Formant update.  Die bug there: fc should be updated, not va.
    // The formants are frozen on a pause phone unless both voice and
    // noise volumes are zero.
    if (tick_208 && (!votraxsc01_locals.rom_pause || !(votraxsc01_locals.filt_fa || votraxsc01_locals.filt_va))) {
        //      interpolate(&votraxsc01_locals.cur_va,  votraxsc01_locals.rom_va);
        interpolate(&votraxsc01_locals.cur_fc, votraxsc01_locals.rom_fc);
        interpolate(&votraxsc01_locals.cur_f1, votraxsc01_locals.rom_f1);
        interpolate(&votraxsc01_locals.cur_f2, votraxsc01_locals.rom_f2);
        interpolate(&votraxsc01_locals.cur_f2q, votraxsc01_locals.rom_f2q);
        interpolate(&votraxsc01_locals.cur_f3, votraxsc01_locals.rom_f3);

        printf("formant int fa=%x va=%x fc=%x f1=%x f2=%02x f2q=%02x f3=%x\n",
               votraxsc01_locals.cur_fa >> 4, votraxsc01_locals.cur_va >> 4,
               votraxsc01_locals.cur_fc >> 4, votraxsc01_locals.cur_f1 >> 4,
               votraxsc01_locals.cur_f2 >> 3, votraxsc01_locals.cur_f2q >> 4,
               votraxsc01_locals.cur_f3 >> 4);
    }

    // Non-formant update. Same bug there, va should be updated, not fc.
    if(tick_625) {

        if (votraxsc01_locals.ticks >= votraxsc01_locals.rom_vd) {
            interpolate(&votraxsc01_locals.cur_fa, votraxsc01_locals.rom_fa);
        }

        if (votraxsc01_locals.ticks >= votraxsc01_locals.rom_cld) {
            //          interpolate(&votraxsc01_locals.cur_fc, votraxsc01_locals.rom_fc);
            interpolate(&votraxsc01_locals.cur_va, votraxsc01_locals.rom_va);
        }

        printf("non-formant int fa=%x va=%x fc=%x f1=%x f2=%02x f2q=%02x f3=%x\n",
               votraxsc01_locals.cur_fa >> 4, votraxsc01_locals.cur_va >> 4, votraxsc01_locals.cur_fc >> 4,
               votraxsc01_locals.cur_f1 >> 4, votraxsc01_locals.cur_f2 >> 3, votraxsc01_locals.cur_f2q >> 4,
               votraxsc01_locals.cur_f3 >> 4);
    }

    // Closure counter, reset every other tick in theory when not
    // active (on the extra rom cycle).
    //
    // The closure level is immediatly used in the analog path,
    // there's no pitch synchronization.

    if (!votraxsc01_locals.cur_closure && (votraxsc01_locals.filt_fa || votraxsc01_locals.filt_va)) {
        votraxsc01_locals.closure = 0;
    }
    else if (votraxsc01_locals.closure != (7 << 2)) {
        votraxsc01_locals.closure++;
    }

    // Pitch counter.  Equality comparison, so it's possible to make
    // it miss by manipulating the inflection inputs, but it'll wrap.
    // There's a delay, hence the +2.

    // Intrinsically pre-divides by two, so we added one bit on the 7

    votraxsc01_locals.pitch = (votraxsc01_locals.pitch + 1) & 0xff;
    if(votraxsc01_locals.pitch == (0xe0 ^ (votraxsc01_locals.inflection << 5) ^ (votraxsc01_locals.filt_f1 << 1)) + 2) {
        votraxsc01_locals.pitch = 0;
    }

    // Filters are updated in index 1 of the pitch wave, which does
    // indeed mean four times in a row.
    if ((votraxsc01_locals.pitch & 0xf9) == 0x08) {
        filters_commit(0);
    }

    // Noise shift register.  15 bits, with a nxor on the last two
    // bits for the loop.
    int inp = (1 || votraxsc01_locals.filt_fa) && votraxsc01_locals.cur_noise && (votraxsc01_locals.noise != 0x7fff);
    votraxsc01_locals.noise = ((votraxsc01_locals.noise << 1) & 0x7ffe) | inp;
    votraxsc01_locals.cur_noise = !(((votraxsc01_locals.noise >> 14) ^ (votraxsc01_locals.noise >> 13)) & 1);

    printf("tick %02x.%03x 625=%d 208=%d pitch=%02x.%x ns=%04x ni=%d noise=%d cl=%x.%x clf=%d/%d ar=%d\n",
           votraxsc01_locals.ticks, votraxsc01_locals.phonetick, tick_625, tick_208, votraxsc01_locals.pitch >> 3,
           votraxsc01_locals.pitch & 7, votraxsc01_locals.noise, inp, votraxsc01_locals.cur_noise,
           votraxsc01_locals.closure >> 2, votraxsc01_locals.closure & 3, votraxsc01_locals.rom_closure,
           votraxsc01_locals.cur_closure, votraxsc01_locals.ar_state);
}

// Shift a history of values by one and insert the new value at the front
static void shift_hist(const double val, double * const hist_array, const size_t N) {
    for(size_t i=N-1; i>0; i--) {
        hist_array[i] = hist_array[i - 1];
    }
    hist_array[0] = val;
}

// Apply a filter and compute the result. 'a' is applied to x (inputs) and 'b' to y (outputs)
static double apply_filter(const double* const x, const double* const y, const double* const a,
                           const size_t Na, const double* const b, const size_t Nb)
{
    double total = 0;
    for(size_t i=0; i<Na; i++) {
        total += x[i] * a[i];
    }
    for(size_t i=1; i<Nb; i++) {
        total -= y[i - 1] * b[i];
    }
    return total / b[0];
}

static float analog_calc()
{
    // Voice-only path.
    // 1. Pick up the pitch wave

    double v = votraxsc01_locals.pitch >= (9 << 3) ? 0 : s_glottal_wave[votraxsc01_locals.pitch >> 3];

    // 2. Multiply by the initial amplifier.  It's linear on the die,
    // even if it's not in the patent.
    v = v * votraxsc01_locals.filt_va * (1.0/15.0);
    shift_hist(v, votraxsc01_locals.voice_1, 4);

    // 3. Apply the f1 filter
    v = apply_filter(votraxsc01_locals.voice_1, votraxsc01_locals.voice_2, votraxsc01_locals.f1_a, 4, votraxsc01_locals.f1_b, 4);
    shift_hist(v, votraxsc01_locals.voice_2, 4);

    // 4. Apply the f2 filter, voice half
    v = apply_filter(votraxsc01_locals.voice_2, votraxsc01_locals.voice_3, votraxsc01_locals.f2v_a, 4, votraxsc01_locals.f2v_b, 4);
    shift_hist(v, votraxsc01_locals.voice_3, 4);

    // Noise-only path
    // 5. Pick up the noise pitch.  Amplitude is linear.  Base
    // intensity should be checked w.r.t the voice.
    double n = 1e4 * ((votraxsc01_locals.pitch & 0x40 ? votraxsc01_locals.cur_noise : 0) ? 1 : -1);
    n = n * votraxsc01_locals.filt_fa * (1.0/15.0);
    shift_hist(n, votraxsc01_locals.noise_1, 3);

    // 6. Apply the noise shaper
    n = apply_filter(votraxsc01_locals.noise_1, votraxsc01_locals.noise_2, votraxsc01_locals.fn_a, 3, votraxsc01_locals.fn_b, 3);
    shift_hist(n, votraxsc01_locals.noise_2, 3);

    // 7. Scale with the f2 noise input
    double n2 = n * votraxsc01_locals.filt_fc * (1.0/15.0);
    shift_hist(n2, votraxsc01_locals.noise_3, 2);

    // 8. Apply the f2 filter, noise half,
    n2 = apply_filter(votraxsc01_locals.noise_3, votraxsc01_locals.noise_4, votraxsc01_locals.f2n_a, 2, votraxsc01_locals.f2n_b, 2);
    shift_hist(n2, votraxsc01_locals.noise_4, 2);

    // Mixed path
    // 9. Add the f2 voice and f2 noise outputs
    double vn = v + n2;
    shift_hist(vn, votraxsc01_locals.vn_1, 4);

    // 10. Apply the f3 filter
    vn = apply_filter(votraxsc01_locals.vn_1, votraxsc01_locals.vn_2, votraxsc01_locals.f3_a, 4, votraxsc01_locals.f3_b, 4);
    shift_hist(vn, votraxsc01_locals.vn_2, 4);

    // 11. Second noise insertion
    vn += n * (5 + (15 ^ votraxsc01_locals.filt_fc)) * (1.0/20.0);
    shift_hist(vn, votraxsc01_locals.vn_3, 4);

    // 12. Apply the f4 filter
    vn = apply_filter(votraxsc01_locals.vn_3, votraxsc01_locals.vn_4, votraxsc01_locals.f4_a, 4, votraxsc01_locals.f4_b, 4);
    shift_hist(vn, votraxsc01_locals.vn_4, 4);

    // 13. Apply the glottal closure amplitude, also linear
    vn = vn * (7 ^ (votraxsc01_locals.closure >> 2)) * (1.0/7.0);
    shift_hist(vn, votraxsc01_locals.vn_5, 2);

    // 13. Apply the final fixed filter
    vn = apply_filter(votraxsc01_locals.vn_5, votraxsc01_locals.vn_6, votraxsc01_locals.fx_a, 1, votraxsc01_locals.fx_b, 2);
    shift_hist(vn, votraxsc01_locals.vn_6, 2);

    return (float)vn * (float)(1.0/2.4); // prevent excessive clipping //!! MAME had a similar magic of * 1.5 here, which is way too loud though
}


/*--------------------------------------------------------------------
 *
 ---------------------------------------------------------------------*/
void Votrax_Update(int16_t* const buffer, const size_t length)
{
    float* const __restrict buffer_f = (float*)buffer;

    printf("Votrax SC-01: update %zu\n", length);

    for (int i = 0; i<length; i++) {
        votraxsc01_locals.sample_count++;
        if (votraxsc01_locals.sample_count & 1) {
            chip_update();
        }
        buffer_f[i] = analog_calc();
        //printf("Votrax SC-01: buffer %d\n", ac);
    }
}



static double global_offset;
static mame_timer* callback_timer = 0;
static double callback_timer_expire_time = 0;
static int callback_timer_modified;

/* list of active timers */
#define MAX_TIMERS 256
static mame_timer timers[MAX_TIMERS];
static mame_timer *timer_head = 0;
static mame_timer *timer_free_head = 0;
static mame_timer *timer_free_tail = 0;

/*-------------------------------------------------
        timer_init - initialize the timer system
-------------------------------------------------*/
void timer_init()
{
    /* we need to wait until the first call to timer_cyclestorun before using real CPU times */
    global_offset = 0.0;
    callback_timer = NULL;
    callback_timer_modified = 0;

    /* reset the timers */
    memset(timers, 0, sizeof(timers));

    /* initialize the lists */
    timer_head = NULL;
    timer_free_head = &timers[0];
    for (int i = 0; i < MAX_TIMERS-1; i++)
    {
        timers[i].tag = -1;
        timers[i].next = &timers[i+1];
    }
    timers[MAX_TIMERS-1].next = NULL;
    timer_free_tail = &timers[MAX_TIMERS-1];
}




/*-------------------------------------------------
        get_relative_time - return the current time
        relative to the global_offset
-------------------------------------------------*/
double get_relative_time()
{
    /* if we're executing as a particular CPU, use its local time as a base */
    /*activecpu = cpu_getactivecpu();
    if (activecpu >= 0)
        return cpunum_get_localtime(activecpu);*/

    /* if we're currently in a callback, use the timer's expiration time as a base */
    if (callback_timer)
        return callback_timer_expire_time;

    /* otherwise, return 0 */
    return 0;
}

/*-------------------------------------------------
        timer_new - allocate a new timer
-------------------------------------------------*/
mame_timer *timer_new()
{
    mame_timer* timer;

    /* remove an empty entry */
    if (!timer_free_head) return NULL;
    timer = timer_free_head;
    timer_free_head = timer->next;
    if (!timer_free_head)
        timer_free_tail = NULL;

    return timer;
}

/*-------------------------------------------------
        timer_list_insert - insert a new timer into
        the list at the appropriate location
-------------------------------------------------*/

void timer_list_insert(mame_timer *timer)
{
    double expire = timer->enabled ? timer->expire : TIME_NEVER;
    mame_timer *t, *lt = NULL;

    /* loop over the timer list */
    for (t = timer_head; t; lt = t, t = t->next)
    {
        /* if the current list entry expires after us, we should be inserted before it */
        /* note that due to floating point rounding, we need to allow a bit of slop here */
        /* because two equal entries -- within rounding precision -- need to sort in */
        /* the order they were inserted into the list */
        if ((t->expire - expire) > TIME_IN_NSEC(1))
        {
            /* link the new guy in before the current list entry */
            timer->prev = t->prev;
            timer->next = t;

            if (t->prev)
                t->prev->next = timer;
            else
                timer_head = timer;
            t->prev = timer;
            return;
        }
    }

    /* need to insert after the last one */
    if (lt)
        lt->next = timer;
    else
        timer_head = timer;
    timer->prev = lt;
    timer->next = NULL;
}

/*-------------------------------------------------
timer_alloc - allocate a permament timer that
        isn't primed yet
-------------------------------------------------*/
mame_timer* timer_alloc(void (*callback)(int))
{
    double time = get_relative_time();
    mame_timer* timer = timer_new();

    /* fail if we can't allocate a new entry */
    if (!timer)
        return NULL;

    /* fill in the record */
    timer->callback = callback;
    timer->callback_param = 0;
    timer->enabled = 0;
    timer->temporary = 0;
    //timer->tag = get_resource_tag();
    timer->period = 0;

    /* compute the time of the next firing and insert into the list */
    timer->start = time;
    timer->expire = TIME_NEVER;
    timer_list_insert(timer);

    /* return a handle */
    return timer;
}
