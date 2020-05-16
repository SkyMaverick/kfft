#define KFFT_ALGO_PLAN_PREPARE(st, flags, type, memneed, pool, len_value)                          \
    do {                                                                                           \
        kfft_pool_t* mm_##type = NULL;                                                             \
        if (lenmem == NULL) {                                                                      \
            mm_##type = (pool) ? pool : kfft_allocator_create(memneed);                            \
            if (mm_##type)                                                                         \
                st = kfft_internal_alloc(mm_##type, sizeof(type));                                 \
            if ((st == NULL) && (pool == NULL))                                                    \
                kfft_allocator_free(mm_##type);                                                    \
        } else {                                                                                   \
            if (pool && *len_value >= (memneed)) {                                                 \
                mm_##type = pool;                                                                  \
                if ((flags)&KFFT_FLAG_RENEW)                                                       \
                    kfft_allocator_clear(mm_##type);                                               \
                st = kfft_internal_alloc(mm_##type, sizeof(type));                                 \
            }                                                                                      \
            *len_value = (memneed);                                                                \
        }                                                                                          \
        if (st)                                                                                    \
            st->object.mmgr = mm_##type;                                                           \
    } while (0)

#define KFFT_ALGO_PLAN_TERMINATE(st, pool)                                                         \
    do {                                                                                           \
        if ((pool) == NULL)                                                                        \
            kfft_allocator_free(KFFT_PLAN_MMGR(st));                                               \
    } while (0)
