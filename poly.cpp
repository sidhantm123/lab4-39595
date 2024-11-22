#include "poly.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>
using namespace std;
using PCPair = std::pair<size_t, int>;

/** Helper function to merge terms with the same power */
void simplify_polynomial(std::vector<PCPair>& terms) {
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

    // Aggregate coefficients
    for (const auto& term : terms) {
        int power = term.first;
        double coeff = term.second;
        aggregated[power] += coeff;
    }

    // Convert the map back to a vector
    std::vector<PCPair> result(aggregated.begin(), aggregated.end());
    terms = result;

}

/** Default constructor: represents the zero polynomial */
polynomial::polynomial() {
    // Add a polynomial of 0x^0 into the group of polynomials that we have
    terms.push_back({0,0});

}

/** Constructor from iterator range */
// template <typename Iter>
// polynomial::polynomial(Iter begin, Iter end) : terms(begin, end) {
//     simplify_polynomial(terms);
// }

/** Copy constructor */
polynomial::polynomial(const polynomial& other) : terms(other.terms) {}

/** Print the polynomial (for debugging) */
void polynomial::print() const {
    for (const auto& [p, c] : terms) {
        std::cout << c << "x^" << p << " ";
    }
    std::cout << std::endl;
}

/** Assignment operator (takes a deep copy of the polynom it is being assigned to)*/
polynomial& polynomial::operator=(const polynomial& other) {
    if (this != &other) {
        terms = other.terms;
    }
    return *this;
}

/** Addition of two polynomials */
polynomial operator+(const polynomial& lhs, const polynomial& rhs) {
    // Get the terms of both polynomials
    std::vector<std::pair<power, coeff>> terms1 = lhs.canonical_form();
    std::vector<std::pair<power, coeff>> terms2 = rhs.canonical_form();

    // Resultant terms
    std::vector<std::pair<power, coeff>> result;

    // Merging process
    size_t i = 0, j = 0;
    while (i < terms1.size() && j < terms2.size()) {
        // if terms have the same power then add coeffs
        if (terms1[i].first == terms2[j].first) {
            // Same degree: add coefficients
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

    // Add any remaining terms from terms1
    while (i < terms1.size()) {
        result.push_back(terms1[i]);
        ++i;
    }

    // Add any remaining terms from terms2
    while (j < terms2.size()) {
        result.push_back(terms2[j]);
        ++j;
    }

    // Simplify the result
    simplify_polynomial(result);

    // Construct a new polynomial with the resultant terms
    return polynomial(result.begin(), result.end());
}

/** Addition of polynomial and integer (poly + int) */
polynomial operator+(const polynomial& poly, int value) {
    // Get the canonical form of the input polynomial
    std::vector<std::pair<power, coeff>> terms = poly.canonical_form();

    // Check if there is a constant term (power = 0)
    bool constantTermFound = false;
    for (auto& [power, coeff] : terms) {
        if (power == 0) {
            coeff += value; // Add the integer to the constant term
            constantTermFound = true;
            break;
        }
    }

    // If no constant term exists, add one with the given value
    // Adding the integer to the end of the poly
    if (!constantTermFound && value != 0) {
        terms.emplace_back(0, value);
    }

    // Sort terms in asc order
    std::sort(terms.begin(), terms.end());

    // Construct and return the resulting polynomial
    return polynomial(terms.begin(), terms.end());
}

/** Addition of integer and polynomial (int + poly) */
polynomial operator+(int value, const polynomial& poly) {
    return poly + value;
}

/** Multiplication of two polynomials */
polynomial operator*(const polynomial& lhs, const polynomial& rhs) {
    std::vector<std::pair<power, coeff>> resultTerms;
    // Get canonical forms of both polynomials
    const auto& lhsTerms = lhs.canonical_form();
    const auto& rhsTerms = rhs.canonical_form();

    // Use a map to accumulate terms by power
    std::map<power, coeff> termMap;

    // Multiply each term in lhs by each term in rhs
    for (const auto& [p1, c1] : lhsTerms) {
        for (const auto& [p2, c2] : rhsTerms) {
            termMap[p1 + p2] += c1 * c2;
        }
    }

    // Convert map back to vector
    for (const auto& [power, coefficient] : termMap) {
        if (coefficient != 0) { // Skip zero coefficients
            resultTerms.emplace_back(power, coefficient);
        }
    }
    // Sort terms by descending power
    std::sort(resultTerms.rbegin(), resultTerms.rend());

    // Construct the result polynomial
    return polynomial(resultTerms.begin(), resultTerms.end());
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
    // Validate divisor is not a zero polynomial
    if (d.canonical_form().empty() || 
        (d.canonical_form().size() == 1 && d.canonical_form()[0].second == 0)) {
        throw std::invalid_argument("Division by zero polynomial");
    }

    // Get canonical forms of the dividend and divisor
    auto dividendTerms = p.canonical_form();
    auto divisorTerms = d.canonical_form();

    // Initialize the dividend and divisor as polynomials
    polynomial dividend(dividendTerms.begin(), dividendTerms.end());
    polynomial divisor(divisorTerms.begin(), divisorTerms.end());

    // Result polynomial to store the remainder
    polynomial remainder = dividend;

    // Perform polynomial division
    while (!remainder.canonical_form().empty() && 
           remainder.canonical_form()[0].first >= divisor.canonical_form()[0].first) {
        // Get the leading terms of the dividend and divisor
        power diff_power = remainder.canonical_form()[0].first - divisor.canonical_form()[0].first;
        coeff quotient_coeff = remainder.canonical_form()[0].second / divisor.canonical_form()[0].second;

        // Create a temporary polynomial for the scaled divisor
        std::vector<std::pair<power, coeff>> tempTerms;
        for (const auto& [p, c] : divisor.canonical_form()) {
            tempTerms.emplace_back(p + diff_power, c * quotient_coeff);
        }
        polynomial temp(tempTerms.begin(), tempTerms.end());

        // Subtract the scaled divisor from the dividend
        remainder = remainder + (temp * -1);
    }

    // Simplify the remainder polynomial
    return remainder;
}

/** Returns the degree of the polynomial */
size_t polynomial::find_degree_of() {
    return terms.empty() ? 0 : terms.rbegin()->first;
}

/** Returns the canonical form of the polynomial */
std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    std::vector<std::pair<power, coeff>> canonical;
    for (auto it = terms.begin(); it != terms.end(); ++it) {
        canonical.emplace_back(it->first, it->second);
    }
    if (canonical.empty()) {
        canonical.emplace_back(0, 0);
    }
    return canonical;
}