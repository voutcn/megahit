#include <pthread.h>
#include <stdlib.h>
#include <limits.h>

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
    struct kt_for_t *t;
    long i;
} ktf_worker_t;

typedef struct kt_for_t {
    int n_threads;
    long n, i;
    ktf_worker_t *w;
    void (*func)(void *, long, int);
    void *data;
} kt_for_t;

static inline long steal_work(kt_for_t *t) {
    int i, min_i = -1;
    long k, min = LONG_MAX;

    for (i = 0; i < t->n_threads; ++i)
        if (min > t->w[i].i) min = t->w[i].i, min_i = i;

    if (t->w[min_i].i < t->n_threads) return -2;

    k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
    return k >= t->n ? -1 : k;
}

static void *ktf_worker(void *data) {
    ktf_worker_t *w = (ktf_worker_t *)data;
    long i;

    for (;;) {
        i = __sync_fetch_and_add(&w->i, w->t->n_threads);

        if (i >= w->t->n) break;

        w->t->func(w->t->data, i, w - w->t->w);
    }

    while ((i = steal_work(w->t)) != -1)
        if (i >= 0) w->t->func(w->t->data, i, w - w->t->w);

    pthread_exit(0);
}

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n) {
    int i;
    kt_for_t t;
    pthread_t *tid;
    t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
    t.w = (ktf_worker_t *)alloca(n_threads * sizeof(ktf_worker_t));
    tid = (pthread_t *)alloca(n_threads * sizeof(pthread_t));

    for (i = 0; i < n_threads; ++i)
        t.w[i].t = &t, t.w[i].i = i;

    for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);

    for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}

static void *ktdf_worker(void *data) {
    ktf_worker_t *w = (ktf_worker_t *)data;
    long i;

    if (w->i < w->t->n) {
        w->t->func(w->t->data, w->i, w->i);
    }

    for (;;) {
        i = __sync_fetch_and_add(&w->t->i, 1);

        if (i >= w->t->n) break;

        w->t->func(w->t->data, i, w->i);
    }

    pthread_exit(0);
}

void kt_dfor(int n_threads, void (*func)(void *, long, int), void *data, long n) {
    int i;
    kt_for_t t;
    pthread_t *tid;
    t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
    t.w = (ktf_worker_t *)alloca(n_threads * sizeof(ktf_worker_t));
    tid = (pthread_t *)alloca(n_threads * sizeof(pthread_t));

    for (i = 0; i < n_threads; ++i)
        t.w[i].t = &t, t.w[i].i = i;

    t.i = n_threads;

    for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktdf_worker, &t.w[i]);

    for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}
