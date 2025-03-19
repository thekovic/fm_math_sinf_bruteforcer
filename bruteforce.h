#ifndef __BRUTEFORCE_H__
#define __BRUTEFORCE_H__

#if BRUTEFORCE == 1
#define BRUTEFORCE_LOOP(id, step, lower_bound, upper_bound, body) \
    int lower_bound_##id = lower_bound; \
    for (int bruteforce_adjust_##id = lower_bound_##id; bruteforce_adjust_##id <= upper_bound; bruteforce_adjust_##id += step) { \
        int32_t as_int_##id = BITCAST_F2I((exact_coefs[id])); \
        int bruteforced_##id = as_int_##id + bruteforce_adjust_##id; \
        float as_float_##id = BITCAST_I2F(bruteforced_##id); \
        bruteforced_coefs[id] = as_float_##id; \
        body \
    }
#else
#define BRUTEFORCE_LOOP(id, step, lower_bound, upper_bound, body) \
{ \
    int lower_bound_##id = lower_bound; \
    int bruteforce_adjust_##id = lower_bound_##id; \
    body \
}
#endif

#endif
