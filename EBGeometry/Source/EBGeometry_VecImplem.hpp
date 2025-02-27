/* EBGeometry
 * Copyright © 2022 Robert Marskar
 * Please refer to Copyright.txt and LICENSE in the EBGeometry root directory.
 */

/*!
  @file   EBGeometry_VecImplem.hpp
  @brief  Implementation of EBGeometry_Vec.hpp
  @author Robert Marskar
*/

#ifndef EBGeometry_VecImplem
#define EBGeometry_VecImplem

// Std includes
#include <math.h>
#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>

// Our includes
#include "EBGeometry_Vec.hpp"
#include "EBGeometry_NamespaceHeader.hpp"

template <typename T>
inline Vec2T<T>::Vec2T()
{
  *this = Vec2T<T>::zero();
}

template <typename T>
inline Vec2T<T>::Vec2T(const Vec2T& u) noexcept
{
  x = u.x;
  y = u.y;
}

template <typename T>
inline constexpr Vec2T<T>::Vec2T(const T& a_x, const T& a_y) : x(a_x), y(a_y)
{}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::zero() noexcept
{
  return Vec2T<T>(T(0.0), T(0.0));
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::one() noexcept
{
  return Vec2T<T>(T(1.0), T(1.0));
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::min() noexcept
{
  return Vec2T<T>(-std::numeric_limits<T>::max(), -std::numeric_limits<T>::max());
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::max() noexcept
{
  return Vec2T<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::infinity() noexcept
{
  return Vec2T<T>(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
}

template <typename T>
inline constexpr Vec2T<T>&
Vec2T<T>::operator=(const Vec2T<T>& u) noexcept
{
  x = u.x;
  y = u.y;

  return (*this);
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::operator+(const Vec2T<T>& u) const noexcept
{
  return Vec2T<T>(x + u.x, y + u.y);
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::operator-(const Vec2T<T>& u) const noexcept
{
  return Vec2T<T>(x - u.x, y - u.y);
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::operator-() const noexcept
{
  return Vec2T<T>(-x, -y);
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::operator*(const T& s) const noexcept
{
  return Vec2T<T>(x * s, y * s);
}

template <typename T>
inline constexpr Vec2T<T>
Vec2T<T>::operator/(const T& s) const noexcept
{
  const T is = 1. / s;

  return Vec2T<T>(x * is, y * is);
}

template <typename T>
inline constexpr Vec2T<T>&
Vec2T<T>::operator+=(const Vec2T<T>& u) noexcept
{
  x += u.x;
  y += u.y;

  return (*this);
}

template <typename T>
inline constexpr Vec2T<T>&
Vec2T<T>::operator-=(const Vec2T<T>& u) noexcept
{
  x -= u.x;
  y -= u.y;

  return (*this);
}

template <typename T>
inline constexpr Vec2T<T>&
Vec2T<T>::operator*=(const T& s) noexcept
{
  x *= s;
  y *= s;

  return (*this);
}

template <typename T>
inline constexpr Vec2T<T>&
Vec2T<T>::operator/=(const T& s) noexcept
{
  const T is = 1. / s;

  x *= is;
  y *= is;

  return (*this);
}

template <typename T>
inline constexpr T
Vec2T<T>::dot(const Vec2T<T>& u) const noexcept
{
  return x * u.x + y * u.y;
}

template <typename T>
inline constexpr T
Vec2T<T>::length() const noexcept
{
  return sqrt(x * x + y * y);
}

template <typename T>
inline constexpr T
Vec2T<T>::length2() const noexcept
{
  return x * x + y * y;
}

template <typename T>
inline constexpr Vec2T<T>
operator*(const T& s, const Vec2T<T>& a_other) noexcept
{
  return a_other * s;
}

template <typename T>
inline constexpr Vec2T<T>
operator/(const T& s, const Vec2T<T>& a_other) noexcept
{
  return a_other / s;
}

template <typename T>
inline Vec3T<T>::Vec3T() noexcept
{
  (*this) = Vec3T<T>::zero();
}

template <typename T>
inline Vec3T<T>::Vec3T(const Vec3T<T>& u) noexcept
{
  m_X[0] = u[0];
  m_X[1] = u[1];
  m_X[2] = u[2];
}

template <typename T>
inline constexpr Vec3T<T>::Vec3T(const T& a_x, const T& a_y, const T& a_z) noexcept : m_X{a_x, a_y, a_z}
{}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::zero() noexcept
{
  return Vec3T<T>(0, 0, 0);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::one() noexcept
{
  return Vec3T<T>(1, 1, 1);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::unit(const size_t a_dir) noexcept
{
  Vec3T<T> v        = Vec3T<T>::zero();
  (v.m_X).at(a_dir) = 1.0;

  return v;
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::min() noexcept
{
  return Vec3T<T>(-std::numeric_limits<T>::max(), -std::numeric_limits<T>::max(), -std::numeric_limits<T>::max());
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::max() noexcept
{
  return Vec3T<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::infinity() noexcept
{
  return Vec3T<T>(
    std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
}

template <typename T>
inline constexpr bool
Vec3T<T>::lessLX(const Vec3T<T>& u) const noexcept
{
  const auto& myComps = std::tuple_cat(m_X);
  const auto& uComps  = std::tuple_cat(u.m_X);

  return std::tie(myComps) < std::tie(uComps);
}

template <typename T>
inline constexpr Vec3T<T>&
Vec3T<T>::operator=(const Vec3T<T>& u) noexcept
{
  m_X[0] = u[0];
  m_X[1] = u[1];
  m_X[2] = u[2];

  return (*this);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator+(const Vec3T<T>& u) const noexcept
{
  return Vec3T<T>(m_X[0] + u[0], m_X[1] + u[1], m_X[2] + u[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator-(const Vec3T<T>& u) const noexcept
{
  return Vec3T<T>(m_X[0] - u[0], m_X[1] - u[1], m_X[2] - u[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator-() const noexcept
{
  return Vec3T<T>(-m_X[0], -m_X[1], -m_X[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator*(const T& s) const noexcept
{
  return Vec3T<T>(s * m_X[0], s * m_X[1], s * m_X[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator*(const Vec3T<T>& s) const noexcept
{
  return Vec3T<T>(s[0] * m_X[0], s[1] * m_X[1], s[2] * m_X[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator/(const T& s) const noexcept
{
  const T is = 1. / s;
  return Vec3T<T>(is * m_X[0], is * m_X[1], is * m_X[2]);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::operator/(const Vec3T<T>& v) const noexcept
{
  return Vec3T<T>(m_X[0] / v[0], m_X[1] / v[1], m_X[2] / v[2]);
}

template <typename T>
inline constexpr Vec3T<T>&
Vec3T<T>::operator+=(const Vec3T<T>& u) noexcept
{
  m_X[0] += u[0];
  m_X[1] += u[1];
  m_X[2] += u[2];

  return (*this);
}

template <typename T>
inline constexpr Vec3T<T>&
Vec3T<T>::operator-=(const Vec3T<T>& u) noexcept
{
  m_X[0] -= u[0];
  m_X[1] -= u[1];
  m_X[2] -= u[2];

  return (*this);
}

template <typename T>
inline constexpr Vec3T<T>&
Vec3T<T>::operator*=(const T& s) noexcept
{
  m_X[0] *= s;
  m_X[1] *= s;
  m_X[2] *= s;

  return (*this);
}

template <typename T>
inline constexpr Vec3T<T>&
Vec3T<T>::operator/=(const T& s) noexcept
{
  const T is = 1. / s;

  m_X[0] *= is;
  m_X[1] *= is;
  m_X[2] *= is;

  return (*this);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::cross(const Vec3T<T>& u) const noexcept
{
  return Vec3T<T>(m_X[1] * u[2] - m_X[2] * u[1], m_X[2] * u[0] - m_X[0] * u[2], m_X[0] * u[1] - m_X[1] * u[0]);
}

template <typename T>
inline constexpr T&
Vec3T<T>::operator[](size_t i) noexcept
{
  return m_X.at(i);
}

template <typename T>
inline constexpr const T&
Vec3T<T>::operator[](size_t i) const noexcept
{
  return m_X.at(i);
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::min(const Vec3T<T>& u) noexcept
{
  m_X[0] = std::min(m_X[0], u[0]);
  m_X[1] = std::min(m_X[1], u[1]);
  m_X[2] = std::min(m_X[2], u[2]);

  return *this;
}

template <typename T>
inline constexpr Vec3T<T>
Vec3T<T>::max(const Vec3T<T>& u) noexcept
{
  m_X[0] = std::max(m_X[0], u[0]);
  m_X[1] = std::max(m_X[1], u[1]);
  m_X[2] = std::max(m_X[2], u[2]);

  return *this;
}

template <typename T>
inline size_t
Vec3T<T>::minDir(const bool a_doAbs) const noexcept
{
  size_t mDir = 0;

  for (size_t dir = 0; dir < 3; dir++) {
    if (a_doAbs) {
      if (std::abs(m_X[dir]) < std::abs(m_X[mDir])) {
        mDir = dir;
      }
    }
    else {
      if (m_X[dir] < m_X[mDir]) {
        mDir = dir;
      }
    }
  }

  return mDir;
}

template <typename T>
inline size_t
Vec3T<T>::maxDir(const bool a_doAbs) const noexcept
{
  size_t mDir = 0;

  for (size_t dir = 0; dir < 3; dir++) {
    if (a_doAbs) {
      if (std::abs(m_X[dir]) > std::abs(m_X[mDir])) {
        mDir = dir;
      }
    }
    else {
      if (m_X[dir] > m_X[mDir]) {
        mDir = dir;
      }
    }
  }

  return mDir;
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator==(const Vec3T<T>& u) const noexcept
{
  return (m_X[0] == u[0] && m_X[1] == u[1] && m_X[2] == u[2]);
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator!=(const Vec3T<T>& u) const noexcept
{
  return !(*this == u);
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator<(const Vec3T<T>& u) const noexcept
{
  return (m_X[0] < u[0] && m_X[1] < u[1] && m_X[2] < u[2]);
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator>(const Vec3T<T>& u) const noexcept
{
  return (m_X[0] > u[0] && m_X[1] > u[1] && m_X[2] > u[2]);
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator<=(const Vec3T<T>& u) const noexcept
{
  return (m_X[0] <= u[0] && m_X[1] <= u[1] && m_X[2] <= u[2]);
}

template <typename T>
inline constexpr bool
Vec3T<T>::operator>=(const Vec3T<T>& u) const noexcept
{
  return (m_X[0] >= u[0] && m_X[1] >= u[1] && m_X[2] >= u[2]);
}

template <typename T>
inline constexpr T
Vec3T<T>::dot(const Vec3T<T>& u) const noexcept
{
  return m_X[0] * u[0] + m_X[1] * u[1] + m_X[2] * u[2];
}

template <typename T>
inline constexpr T
Vec3T<T>::length() const noexcept
{
  return sqrt(m_X[0] * m_X[0] + m_X[1] * m_X[1] + m_X[2] * m_X[2]);
}

template <typename T>
inline constexpr T
Vec3T<T>::length2() const noexcept
{
  return m_X[0] * m_X[0] + m_X[1] * m_X[1] + m_X[2] * m_X[2];
}

template <typename T>
inline Vec2T<T>
min(const Vec2T<T>& u, const Vec2T<T>& v) noexcept
{
  return Vec2T<T>(std::min(u.x, v.x), std::min(u.y, v.y));
}

template <typename T>
inline Vec2T<T>
max(const Vec2T<T>& u, const Vec2T<T>& v) noexcept
{
  return Vec2T<T>(std::max(u.x, v.x), std::max(u.y, v.y));
}

template <typename T>
inline T
dot(const Vec2T<T>& u, const Vec2T<T>& v) noexcept
{
  return u.dot(v);
}

template <typename T>
inline T
length(const Vec2T<T>& v) noexcept
{
  return v.length();
}

template <class R, typename T>
inline constexpr Vec3T<T>
operator*(const R& s, const Vec3T<T>& a_other) noexcept
{
  return a_other * s;
}

template <typename T>
inline constexpr Vec3T<T>
operator*(const Vec3T<T>& u, const Vec3T<T>& v) noexcept
{
  return u * v;
}

template <class R, typename T>
inline constexpr Vec3T<T>
operator/(const R& s, const Vec3T<T>& a_other) noexcept
{
  return a_other / s;
}

template <typename T>
inline constexpr Vec3T<T>
min(const Vec3T<T>& u, const Vec3T<T>& v) noexcept
{
  return Vec3T<T>(std::min(u[0], v[0]), std::min(u[1], v[1]), std::min(u[2], v[2]));
}

template <typename T>
inline constexpr Vec3T<T>
max(const Vec3T<T>& u, const Vec3T<T>& v) noexcept
{
  return Vec3T<T>(std::max(u[0], v[0]), std::max(u[1], v[1]), std::max(u[2], v[2]));
}

template <typename T>
inline constexpr T
dot(const Vec3T<T>& u, const Vec3T<T>& v) noexcept
{
  return u.dot(v);
}

template <typename T>
inline constexpr T
length(const Vec3T<T>& v) noexcept
{
  return v.length();
}

#include "EBGeometry_NamespaceFooter.hpp"

#endif
