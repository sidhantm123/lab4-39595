#include "poly.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>
using namespace std;
using PCPair = std::pair<size_t, int>;

//use this is condense the polynomial (eliminate 0 terms and combine like terms)
void polynomial::simplify_polynomial(std::vector<PCPair>& terms) {
    int i = 0;
    for (auto temp = terms.begin(); temp != terms.end();) {
        if (temp[i].second == 0) {
            temp = terms.erase(temp);
        } 
        else {
            ++temp;
        }
        i++;
    }
    // create a map to find like terms fast
    std::map<size_t, int> aggregated;

    // add like terms
    for (const auto& term : terms) {
        int power = term.first;
        double coeff = term.second;
        aggregated[power] += coeff;
    }

    // Convert map back to a vector
    std::vector<PCPair> result(aggregated.begin(), aggregated.end());
    terms = result;

}

// default constructor
polynomial::polynomial() {
    terms.push_back({0,0});

}

// iter constructor
// template <typename Iter>
// polynomial::polynomial(Iter begin, Iter end) : terms(begin, end) {
//     simplify_polynomial(terms);
// }

/** Copy constructor */
polynomial::polynomial(const polynomial& other) : terms(other.terms) {}

// Print the polynomial (for debugging)
void polynomial::print() const {
    for (const auto& [p, c] : terms) {
        std::cout << c << "x^" << p << " ";
    }
    std::cout << std::endl;
}

// Assignment operator (takes a deep copy of the polynom it is being assigned to)
polynomial& polynomial::operator=(const polynomial& other) {
    if (this != &other) {
        terms = other.terms;
    }
    return *this;
}

// Polynomial + polynomial
polynomial polynomial::operator+(const polynomial& rhs) {
    // conv to canonical form to traverse and index more easilu
    auto lhs = *this;
    std::vector<PCPair> terms1 = lhs.canonical_form();
    std::vector<PCPair> terms2 = rhs.canonical_form();

    // Result
    std::vector<std::pair<power, coeff>> result;

    // Merge
    size_t i = 0;
    size_t j = 0;
    while (i < terms1.size() && j < terms2.size()) {
        // if terms have the same power then add coeffs
        if (terms1[i].first == terms2[j].first) {
            int coefficientSum = terms1[i].second + terms2[j].second;
            if (coefficientSum != 0) {
                result.emplace_back(terms1[i].first, coefficientSum);
            }
            ++i;
            ++j;
        } else if (terms1[i].first > terms2[j].first) { 
            // if lhs term power is greater it goes in before
            result.push_back(terms1[i]);
            ++i;
        } else {
            // add remaining terms from rhs
            result.push_back(terms2[j]);
            ++j;
        }
    }

    // Add remaining terms from terms1
    while (i < terms1.size()) {
        result.push_back(terms1[i]);
        ++i;
    }

    // Add remaining terms from terms2
    while (j < terms2.size()) {
        result.push_back(terms2[j]);
        ++j;
    }

    polynomial::simplify_polynomial(result);

    // use constructor to convers result to type polynomial and retirn it
    return polynomial(result.begin(), result.end());
}

// poly + int
polynomial polynomial::operator+(const int value) {
    // Get the canonical form
    std::vector<PCPair> terms = this->canonical_form();

    // Check if there is a constant term (power = 0)
    bool constantTermFound = false;
    for (auto& [power, coeff] : terms) {
        if (power == 0) {
            coeff += value;
            constantTermFound = true;
            break;
        }
    }

    // Adding the integer to the end of the poly if there is no other const term
    if (!constantTermFound && value != 0) {
        terms.emplace_back(0, value);
    }
    // Sort terms in asc order
    std::sort(terms.begin(), terms.end());

    return polynomial(terms.begin(), terms.end());
}

// int + poly (just call poly + int)
polynomial operator+(int value, polynomial& poly) {
    return poly + value;
}

// polynomial * polynomial
polynomial polynomial::operator*(const polynomial& rhs) {
    // to store result
    std::vector<PCPair> resultTerms;

    // Get canonical form
    const auto& lhsTerms = this->canonical_form();
    const auto& rhsTerms = rhs.canonical_form();

    // map to traverse and index better
    std::map<power, coeff> termMap;

    // Multiply each term in lhs by each term in rhs
    for (const auto& [p1, c1] : lhsTerms) {
        for (const auto& [p2, c2] : rhsTerms) {
            termMap[p1 + p2] += c1 * c2;
        }
    }

    // Convert map back to vector
    for (const auto& [power, coefficient] : termMap) {
        if (coefficient != 0) {
            resultTerms.emplace_back(power, coefficient);
        }
    }
    // Sort terms asc
    std::sort(resultTerms.begin(), resultTerms.end());

    // Convert std::vector<PCPair> back to polynoomail
    auto toReturn = polynomial(resultTerms.begin(), resultTerms.end());
    return toReturn;
}

// poly * int
polynomial polynomial::operator*(const int value) {
    polynomial result = *this;
    for (auto& [p, c] : result.terms) {
        c *= value;
    }
    simplify_polynomial(result.terms);
    return result;
}

// int * poly (just return poly * int)
polynomial operator*(const int value, polynomial &poly) {
    return poly * value;
}

/** Polynomial modulus (poly % poly) */
polynomial polynomial::operator%(const polynomial& d) {
    const auto p = *this;
    // check if divisor is == 0
    if (d.canonical_form().empty() || 
        (d.canonical_form().size() == 1 && d.canonical_form()[0].second == 0)) {
        throw std::invalid_argument("Division by zero polynomial");
    }

    // Get canonical forms of the dividend and divisor
    auto dividendTerms = p.canonical_form();
    auto divisorTerms = d.canonical_form();

    polynomial dividend = *this;
    polynomial divisor = d;

    // store the remainder
    polynomial remainder = dividend;

    // division
    while (!remainder.canonical_form().empty() && 
           remainder.canonical_form()[0].first >= divisor.canonical_form()[0].first) {
        // Get the leading terms of the dividend and divisor
        int diff_power = remainder.canonical_form()[0].first - divisor.canonical_form()[0].first;
        size_t quotient_coeff = remainder.canonical_form()[0].second / divisor.canonical_form()[0].second;

        // temp vector
        std::vector<PCPair> tempTerms;
        for (const auto& [p, c] : divisor.canonical_form()) {
            tempTerms.emplace_back(p + diff_power, c * quotient_coeff);
        }
        // conv temp to polynomial
        polynomial temp(tempTerms.begin(), tempTerms.end());
        remainder = remainder + (temp * -1);
    }

    polynomial::simplify_polynomial(remainder.terms);
    return remainder;
}

// get degree
size_t polynomial::find_degree_of() {
    return terms.empty() ? 0 : terms.rbegin()->first;
}

// this will give you the canonical form of the polynomial
std::vector<PCPair> polynomial::canonical_form() const {
    std::vector<PCPair> canonical;
    for (auto it = terms.begin(); it != terms.end(); ++it) {
        canonical.emplace_back(it->first, it->second);
    }
    if (canonical.empty()) {
        canonical.emplace_back(0, 0);
    }
    return canonical;
}