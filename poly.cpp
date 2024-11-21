#include "poly.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>

/** Helper function to merge terms with the same power */
void simplify_polynomial(std::map<power, coeff>& terms) {
    for (auto temp = terms.begin(); temp != terms.end();) {
        if (temp->second == 0) {
            temp = terms.erase(temp);
        } else {
            ++temp;
        }
    }
}

/** Default constructor: represents the zero polynomial */
polynomial::polynomial() : terms({{0, 0}}) {}

/** Constructor from iterator range */
template <typename Iter>
polynomial::polynomial(Iter begin, Iter end) : terms(begin, end) {
    simplify_polynomial(terms);
}

/** Copy constructor */
polynomial::polynomial(const polynomial& other) : terms(other.terms) {}

/** Print the polynomial (for debugging) */
void polynomial::print() const {
    for (const auto& [p, c] : terms) {
        std::cout << c << "x^" << p << " ";
    }
    std::cout << std::endl;
}

/** Assignment operator */
polynomial& polynomial::operator=(const polynomial& other) {
    if (this != &other) {
        terms = other.terms;
    }
    return *this;
}

/** Addition of two polynomials */
polynomial operator+(const polynomial& lhs, const polynomial& rhs) {
    polynomial result(lhs);
    for (const auto& [p, c] : rhs.terms) {
        result.terms[p] += c;
    }
    simplify_polynomial(result.terms);
    return result;
}

/** Addition of polynomial and integer (poly + int) */
polynomial operator+(const polynomial& poly, int value) {
    polynomial result(poly);
    result.terms[0] += value;
    simplify_polynomial(result.terms);
    return result;
}

/** Addition of integer and polynomial (int + poly) */
polynomial operator+(int value, const polynomial& poly) {
    return poly + value;
}

/** Multiplication of two polynomials */
polynomial operator*(const polynomial& lhs, const polynomial& rhs) {
    polynomial result;
    result.terms.clear();
    for (const auto& [p1, c1] : lhs.terms) {
        for (const auto& [p2, c2] : rhs.terms) {
            result.terms[p1 + p2] += c1 * c2;
        }
    }
    simplify_polynomial(result.terms);
    return result;
}

/** Multiplication of polynomial and integer (poly * int) */
polynomial operator*(const polynomial& poly, int value) {
    polynomial result(poly);
    for (auto& [p, c] : result.terms) {
        c *= value;
    }
    simplify_polynomial(result.terms);
    return result;
}

/** Multiplication of integer and polynomial (int * poly) */
polynomial operator*(int value, const polynomial& poly) {
    return poly * value;
}

/** Polynomial modulus (poly % poly) */
polynomial operator%(const polynomial& p, const polynomial& d) {
    if (d.terms.empty() || (d.terms.size() == 1 && d.terms.begin()->second == 0)) {
        throw std::invalid_argument("Division by zero polynomial");
    }

    polynomial dividend(p);
    polynomial divisor(d);
    polynomial remainder;

    while (!dividend.terms.empty() && dividend.terms.rbegin()->first >= divisor.terms.rbegin()->first) {
        power diff_power = dividend.terms.rbegin()->first - divisor.terms.rbegin()->first;
        coeff quotient_coeff = dividend.terms.rbegin()->second / divisor.terms.rbegin()->second;
        polynomial temp = divisor * quotient_coeff;
        for (auto& [p, c] : temp.terms) {
            temp.terms[p + diff_power] = c;
        }
        dividend = dividend + temp * -1;
    }
    simplify_polynomial(dividend.terms);
    return dividend;
}

/** Returns the degree of the polynomial */
size_t polynomial::find_degree_of() {
    return terms.empty() ? 0 : terms.rbegin()->first;
}

/** Returns the canonical form of the polynomial */
std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    std::vector<std::pair<power, coeff>> canonical;
    for (auto it = terms.rbegin(); it != terms.rend(); ++it) {
        canonical.emplace_back(it->first, it->second);
    }
    if (canonical.empty()) {
        canonical.emplace_back(0, 0);
    }
    return canonical;
}
