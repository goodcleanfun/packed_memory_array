#include <stdint.h>
#include "greatest/greatest.h"

#define PMA_NAME int_pma
#define PMA_TYPE int32_t
#include "packed_memory_array.h"
#undef PMA_NAME
#undef PMA_TYPE

TEST test_pma(void) {
    int_pma *v = int_pma_new();

    ASSERT(int_pma_insert(v, 5));
    ASSERT(int_pma_insert(v, 3));
    ASSERT(int_pma_insert(v, 8));
    ASSERT(int_pma_insert(v, 1));
    ASSERT(int_pma_insert(v, 7));
    ASSERT(int_pma_insert(v, 2));
    ASSERT(int_pma_delete(v, 5));
    ASSERT(int_pma_insert(v, 6));
    ASSERT(int_pma_insert(v, 4));
    ASSERT(int_pma_delete(v, 6));
    ASSERT(int_pma_insert(v, 9));

    int64_t i, j;
    ASSERT(int_pma_find(v, 1, &i));

    for (int32_t val = 2; val <= 9; val++) {
        if ((val < 5 || val > 6)) {
            ASSERT(int_pma_find(v, val, &j));
            ASSERT(j > i);
            i = j;
        } else {
            ASSERT(!int_pma_find(v, val, &j));
        }
    }
    int_pma_destroy(v);
    PASS();
}

/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();

int32_t main(int32_t argc, char **argv) {
    GREATEST_MAIN_BEGIN();      /* command-line options, initialization. */

    RUN_TEST(test_pma);

    GREATEST_MAIN_END();        /* display results */
}
