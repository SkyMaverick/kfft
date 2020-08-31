/*!
    \file
    \brief Standart predefined algorithms.
 */

/*!
    Plan configure operation basis macro. See ::kfft_config_cpx
    \param[in] plan - plan object (::kfft_plan_cpx for example). Must be NULL.
    \param[in] flags - operation flags ::kfft_eval_flags
    \param[in] type - plan object type (kfft_plan_cpx, kfft_plan_scr, ...)
    \param[in] memneed - needed memory size (bytes)
    \param[in] pool - memory allocator (::kfft_config_cpx)
    \param[in] len_value - user context variable (::kfft_config_cpx)
 */
#define KFFT_ALGO_PLAN_PREPARE(plan, flags, type, memneed, pool, len_value)                        \
    do {                                                                                           \
        kfft_pool_t* mm_##type = NULL;                                                             \
        if (len_value == NULL) {                                                                   \
            mm_##type = (pool) ? pool : kfft_pool_create(memneed);                                 \
            if (mm_##type)                                                                         \
                plan = kfft_pool_alloc(mm_##type, sizeof(type));                                   \
            if ((plan == NULL) && (pool == NULL))                                                  \
                kfft_pool_free(mm_##type);                                                         \
        } else {                                                                                   \
            if (pool && *len_value >= (memneed)) {                                                 \
                mm_##type = pool;                                                                  \
                if ((flags)&KFFT_FLAG_RENEW)                                                       \
                    kfft_pool_clear(mm_##type);                                                    \
                plan = kfft_pool_alloc(mm_##type, *len_value);                                     \
            }                                                                                      \
            if (plan == NULL)                                                                      \
                *len_value = (memneed);                                                            \
        }                                                                                          \
        if (plan)                                                                                  \
            plan->object.mmgr = mm_##type;                                                         \
    } while (0)

/*!
    Terminate and free new allocator if plan don't create.
    \param[in] plan - plan object (::kfft_plan_cpx for example)
    \param[in] pool - memory allocator (::kfft_config_cpx)
 */

#define KFFT_ALGO_PLAN_TERMINATE(plan, pool)                                                       \
    do {                                                                                           \
        if ((pool) == NULL)                                                                        \
            kfft_pool_free(KFFT_PLAN_MMGR(plan));                                                  \
    } while (0)
