#include "fmatha.h"
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#define M_PI (3.14159265358979323846f)

#define MAX_APPROX (5)

static const float pi_hi            = 3.14159274e+00f; // 0x1.921fb6p+01
static const float pi_lo            =-8.74227766e-08f; // -0x1.777a5cp-24
static const float half_pi_hi       =  1.57079637e+0f; //  0x1.921fb6p+0
static const float half_pi_lo       = -4.37113883e-8f; // -0x1.777a5cp-25

#define POLYGON_COEFS (6)

static const float original_coefs[POLYGON_COEFS] = {
    - 1.01321176e-1f,
    6.62087463e-3f,
    - 1.73503853e-4f,
    2.52223435e-6f,
    - 2.33177868e-8f,
    1.32729383e-10f
};

static const float exact_coefs[POLYGON_COEFS] = {
    -0.101321183346709072589001712988183609230944236760490476f,
    0.00662087952180793343258682906697112938547424931632185616f,
    -0.000173505057912483501491115906801116298084629719204655552f,
    2.52229235749396866288379170129828403876289663605034418e-6f,
    -2.33177897192836082466066115718536782354224647348350113e-8f,
    1.32913446369766718120324917415992976452154154051525892e-10f
};

static float bruteforced_coefs[POLYGON_COEFS] = {
    -0.101321183346709072589001712988183609230944236760490476f,
    0.00662087952180793343258682906697112938547424931632185616f,
    -0.000173505057912483501491115906801116298084629719204655552f,
    2.52229235749396866288379170129828403876289663605034418e-6f,
    -2.33177897192836082466066115718536782354224647348350113e-8f,
    1.32913446369766718120324917415992976452154154051525892e-10f
};

__attribute__((noinline))
static float sinf_approx(float x, int approx, float* tested_set) {
    // Approximation of sine to 5 ULP with Chebyshev polynomials
    // http://mooooo.ooo/chebyshev-sine-approximation/
    float p, s;

    p = 0;
    s = x * x;
    // Execute only a portion of the series, depending on the approximation level.
    // This generate the most efficient code among similar approaches.
    if ((--approx < 0)) p += tested_set[5], p *= s;
    if ((--approx < 0)) p += tested_set[4], p *= s;
    if ((--approx < 0)) p += tested_set[3], p *= s;
    if ((--approx < 0)) p += tested_set[2], p *= s;
    if ((--approx < 0)) p += tested_set[1], p *= s;
    if ((--approx < 0)) p += tested_set[0];
    return x * ((x - pi_hi) - pi_lo) * ((x + pi_hi) + pi_lo) * p;   
}

float fm_sinf_approx(float x, int approx, float* tested_set) {
    // sinf_approx has been designed to operate in the [-π, +π] range, so
    // bring the argument there. This reduction using fm_fmodf is not
    // very accurate for large numbers, so it will introduce more error compared
    // to the 5 ULP figure.
    x = fm_fmodf(x+pi_hi, 2*pi_hi) - pi_hi;
    x = sinf_approx(x, approx, tested_set);
    
    return x;
}

float fm_atan2f(float y, float x) {
    // Approximation of atan2f using a polynomial minmax approximation in [0,1]
    // calculated via the Remez algorithm (https://math.stackexchange.com/a/1105038).
    // The reported error is 6.14e-4, so it's precise for at least three decimal
    // digits which is usually more than enough for angles.
    float ay = fabsf(y);
    float ax = fabsf(x);
    float a = (ay < ax) ? ay/ax : ax/ay;
    float s = a * a;
    float r = ((-0.0464964749f * s + 0.15931422f) * s - 0.327622764f) * s * a + a;
    if (ay > ax)
        r = half_pi_hi - r;
    if (BITCAST_F2I(x) < 0) r = pi_hi - r;
    return copysignf(r, y);
}

void print_coefs(float* coefs)
{
    for (int coef = 0; coef < POLYGON_COEFS; coef++)
    {
        int32_t as_int = BITCAST_F2I((coefs[coef]));
        printf("coef %i: %.20f (0x%x) (%i)\n", coef, coefs[coef], as_int, as_int);
    }
}

void run_fm_sinf_over_all_f32s(float* tested_set)
{
    int approx = 0;
    int64_t out_of_bounds_errors = 0;
    double result_diff_sum = 0;
    double max_measured_error = 0;

    for (int64_t as_int = 0; as_int <= 0xffffffff; as_int++)
    {
        float as_float = BITCAST_I2F(as_int);
        if ((as_float > (1 * M_PI)) || (as_float < (-1 * M_PI)))
        {
            continue;
        }
        
        float result = fm_sinf_approx(as_float, approx, tested_set);
        float std_result = sinf(as_float);
        float diff = result - std_result;

        double diff_abs = fabs(diff);
        max_measured_error = (diff_abs > max_measured_error) ? diff_abs : max_measured_error;

        float diff_squared = diff * diff;
        result_diff_sum += diff_squared;

        if (result > 1.0f || result < -1.0f)
        {
            //out_of_bounds_errors++;
            return;
        }
    }

    print_coefs(tested_set);    

    //printf("Out of bounds errors: %li\n", out_of_bounds_errors);
    printf("RMSD: %.20f\n", sqrt(result_diff_sum));
    printf("maximum measured error: %.20f\n", max_measured_error);
    printf("------------------------------------\n");
}

static int current_thread_count = 0;

void* bruteforce_thread_func(void* arg) {
    float* tested_set = (float*) arg;
    run_fm_sinf_over_all_f32s(tested_set);
    free(arg);
    current_thread_count--;
    return NULL;
}

#if 1
#define BRUTEFORCE_LOOP(id, lower_bound, upper_bound, body) \
    for (int bruteforce_adjust_##id = lower_bound; bruteforce_adjust_##id <= upper_bound; bruteforce_adjust_##id++) { \
        int as_int_##id = BITCAST_F2I((exact_coefs[id])); \
        int bruteforced_##id = as_int_##id + bruteforce_adjust_##id; \
        float as_float_##id = BITCAST_I2F(bruteforced_##id); \
        bruteforced_coefs[id] = as_float_##id; \
        body \
    }
#else
#define BRUTEFORCE_LOOP(id, body) { body }
#endif

#define MAX_THREADS (10)

int main()
{
    float* tested_set = bruteforced_coefs;
    printf("Initial coefs:\n");
    print_coefs(tested_set);
    printf("------------------------------------\n");

    int current_thread = 0;
    pthread_t threads[MAX_THREADS * 2];

    for (int approx = 0; approx <= 0; approx++)
    {
        printf("fm_sinf_approx; approx level %i:\n", approx);
        BRUTEFORCE_LOOP(0, -10, 10,
            BRUTEFORCE_LOOP(1, -10, 10,
                BRUTEFORCE_LOOP(2, -10, 10,
                    BRUTEFORCE_LOOP(3, -10, 10,
                        BRUTEFORCE_LOOP(4, -10, 10,
                            BRUTEFORCE_LOOP(5, -10, 10,
                                {
                                    while (current_thread_count > MAX_THREADS)
                                    {
                                        sleep(1);
                                    }

                                    float* current_set = malloc(sizeof(float) * 6);
                                    memcpy(current_set, bruteforced_coefs, sizeof(float) * 6);
                                    printf("%i %i %i\n", bruteforce_adjust_3, bruteforce_adjust_4, bruteforce_adjust_5);
                                    int error_num = pthread_create(&threads[current_thread % MAX_THREADS], NULL, bruteforce_thread_func, current_set);
                                    if (error_num) {
                                        printf("Failed to create thread %i with error code %i\n", current_thread, error_num);
                                        return 1;
                                    }

                                    current_thread_count++;
                                    current_thread++;
                                }
                            )
                        )
                    )
                )
            )
        )
    }

    return 0;
}
