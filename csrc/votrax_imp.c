#include <math.h>
#include <memory.h>
#include <stdio.h>

#include "votrax_imp.h"

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
    int i;

    // Only these two counters are reset on phone change, the rest is
    // free-running.
    votraxsc01_locals.phonetick = 0;
    votraxsc01_locals.ticks = 0;

    // In the real chip, the rom is re-read all the time.  Since it's
    // internal and immutable, no point in not caching it though.
    for (i = 0; i<64; i++) {
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

int votrax_start()
{
    memset(&votraxsc01_locals, 0x00, sizeof(votraxsc01_locals));

    // Defines the sound interface
    //votraxsc01_locals.intf = msound->sound_interface;

    // initialize internal state
    // TODO: Figure out what a good starting clock is
    votraxsc01_locals.mainclock = 9000; //votraxsc01_locals.intf->baseFrequency[0]; //!! clock();
    votraxsc01_locals.sclock = votraxsc01_locals.mainclock / 18.0;
    votraxsc01_locals.cclock = votraxsc01_locals.mainclock / 36.0;

    // This is used to initialize the device
    //votraxsc01_locals.stream = stream_init_float("Votrax - SC01", votraxsc01_locals.intf->mixing_level[0], votraxsc01_locals.sclock, 0, Votrax_Update, 1);
    //votraxsc01_locals.timer = timer_alloc(VOTRAXSC01_sh_start_timeout);

    // reset outputs
    //!! m_ar_cb.resolve_safe();
    //votraxsc01_locals.ar_state = ASSERT_LINE;

    //!! was the separate reset code from here on:

    // Technically, there's no reset in this chip, and initial state
    // is random.  Still, it's a good idea to start it with something
    // sane.

    votraxsc01_locals.phone = 0x3f;
    votraxsc01_locals.inflection = 0;
    //votraxsc01_locals.ar_state = ASSERT_LINE;
    //!! m_ar_cb(votraxsc01_locals.ar_state);

    votraxsc01_locals.sample_count = 0;

    // Initialize the m_rom* values
    phone_commit();

    // Clear the interpolation sram
    votraxsc01_locals.cur_fa = votraxsc01_locals.cur_fc = votraxsc01_locals.cur_va = 0;
    votraxsc01_locals.cur_f1 = votraxsc01_locals.cur_f2 = votraxsc01_locals.cur_f2q = votraxsc01_locals.cur_f3 = 0;

    // Initialize the m_filt* values and the filter coefficients
    //filters_commit(1);

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
