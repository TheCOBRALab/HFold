#ifndef ARRAY
#define ARRAY

#ifndef VRNA_ARRAY_GROW_FORMULA
/**
 *  @brief The default growth formula for array
 */
#define VRNA_ARRAY_GROW_FORMULA(n)                      (1.4 * (n) + 8)
#endif

/**
 *  @brief  The header of an array
 */
typedef struct vrna_array_header_s {
  size_t  num;  /**< @brief The number of elements in an array */
  size_t  size; /**< @brief The actual capacity of an array */
} vrna_array_header_t;

/**
 *  @brief Define an array
 */
#define vrna_array(Type) Type *

/**
 *  @brief  Retrieve a pointer to the header of an array @p input
 */
#define VRNA_ARRAY_HEADER(input)                        ((vrna_array_header_t *)(input) - 1)
/**
 *  @brief  Get the number of elements of an array @p input
 */
#define vrna_array_size(input)                          (VRNA_ARRAY_HEADER(input)->num)

/**
 *  @brief Initialize an array @p a with a particular pre-allocated size @p init_size
 *
 */
#define vrna_array_init_size(a, init_size) do { \
  void **a_ptr = (void **)&(a); \
  size_t size = sizeof(*(a)) * (init_size) + sizeof(vrna_array_header_t); \
  vrna_array_header_t *h = (vrna_array_header_t *)vrna_alloc(size); \
  h->num           = 0; \
  h->size          = init_size; \
  *a_ptr           = (void *)(h + 1); \
} while (0)

/**
 *  @brief Initialize an array @p a
 */
#define vrna_array_init(a)  vrna_array_init_size(a, VRNA_ARRAY_GROW_FORMULA(0));





#endif