/**
 * @file util.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains utility macro.
 * @date 2024-03-05
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __UTIL_UTIL_H__
#define __UTIL_UTIL_H__

/**
 * @brief Minimum of 2 values.
 *
 * @param _x first value.
 * @param _y second value.
 * @return min(_x, _y).
 *
 * @warning Expression with side-effects multiply evaluated by macro.
 */
#define MIN(_x, _y) (((_x) < (_y)) ? (_x) : (_y))

#endif