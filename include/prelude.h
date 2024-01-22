/**
 * @file prelude.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains prelude to the entire project.
 * @date 2024-01-17
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_PRELUDE_H__
#define __REED_SOLOMON_PRELUDE_H__

/**
 * @brief Reed-Solomon code length: 2^16 - 1.
 */
#define N 65535

/**
 * @brief Symbol (message) size.
 */
#define SYMBOL_SIZE 1300

#ifndef NDEBUG
#define _static
#else
#define _static static
#endif

#endif