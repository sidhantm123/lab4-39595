#ifndef POLY_H
#define POLY_H

#include <vector>
#include <utility>
#include <map>
#include <cstddef>
#include <stdexcept>

using power = size_t;
using coeff = int;

class polynomial
{
private:
    // Internal representation of the polynomial as a map of power to coefficient
    std::map<power, coeff> terms;

public:
    /**
     * @brief Construct a new polynomial object that is the number 0 (i.e., 0x^0)
     */
    polynomial();

    /**
     * @brief Construct a new polynomial object from an iterator to pairs of <power, coeff>
     *
     * @tparam Iter
     *  An iterator that points to a std::pair<power, coeff>
     * @param begin
     *  The start of the container to copy elements from
     * @param end
     *  The end of the container to copy elements from
     */
    template <typename Iter>
    polynomial(Iter begin, Iter end) : terms(begin, end) {
        simplify_polynomial(terms);
    }

    /**
     * @brief Construct a new polynomial object from an existing polynomial object
     *
     * @param other
     *  The polynomial to copy
     */
    polynomial(const polynomial &other);

    /**
     * @brief Prints the polynomial (used for debugging, not graded)
     */
    void print() const;

    /**
     * @brief Assignment operator to copy another polynomial
     *
     * @param other
     * The polynomial to copy
     * @return
     * A reference to the copied polynomial
     */
    polynomial &operator=(const polynomial &other);

    /**
     * @brief Returns the degree of the polynomial
     *
     * @return size_t
     *  The degree of the polynomial
     */
    size_t find_degree_of();

    /**
     * @brief Returns a vector that contains the polynomial in canonical form.
     *        This means that the power at index 0 is the largest power in the polynomial,
     *        the power at index 1 is the second largest power, etc.
     *
     *        Example:
     *        x^2 + 7x^4 + 1 would be returned as [(4, 7), (2, 1), (0, 1)]
     *
     *        Note: Any terms with a coefficient of zero are not returned in the
     *        canonical form, except for the polynomial 0.
     *        For example:
     *        y = 0 would be returned as [(0, 0)]
     *
     * @return std::vector<std::pair<power, coeff>>
     *  A vector of pairs representing the canonical form of the polynomial
     */
    std::vector<std::pair<power, coeff>> canonical_form() const;

    /**
     * @brief Overload the +, *, and % operators for polynomial arithmetic
     */
    friend polynomial operator+(const polynomial &lhs, const polynomial &rhs);
    friend polynomial operator+(const polynomial &poly, int value);
    friend polynomial operator+(int value, const polynomial &poly);

    friend polynomial operator*(const polynomial &lhs, const polynomial &rhs);
    friend polynomial operator*(const polynomial &poly, int value);
    friend polynomial operator*(int value, const polynomial &poly);

    friend polynomial operator%(const polynomial &p, const polynomial &d);

private:
    /**
     * @brief Simplifies the polynomial by removing terms with a coefficient of 0
     *
     * @param terms
     * A map of power to coefficient representing the polynomial
     */
    static void simplify_polynomial(std::map<power, coeff>& terms);
};

#endif