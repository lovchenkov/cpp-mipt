#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <array>

class BigInteger {
 private:
  bool sign_;  // true iff BigInteger >= 0
  std::vector<int64_t> digits_;  // keeps digits in the reversed order
  void increase();  // *this supposed to be non-negative;
  void decrease();  // --//--

 public:
  const static int base = 100000;
  const static int log_base = 5;
  BigInteger() : sign_(true), digits_(std::vector<int64_t>(0)) {}
  explicit operator bool() { return length() > 0; }
  BigInteger(int64_t number) : sign_(number >= 0) {
    if (number < 0) number *= -1;
    while (number > 0) {
      digits_.push_back(number % base);
      number /= base;
    }
  }
  BigInteger(const std::string& number) {
    int64_t first_digit = (number[0] == '-') ? 1 : 0;
    if (number[first_digit] == '0') {
      sign_ = true;
      digits_ = std::vector<int64_t>(0);
      return;
    }
    for (int i = number.length() - 1; i >= first_digit; i -= log_base) {
      std::string current_digit = "";
      for (int j = std::max(i - log_base + 1, static_cast<int>(first_digit)); j <= i; ++j) {
        current_digit.push_back(number[j]);
      }
      digits_.push_back(std::stoi(current_digit));
    }
    sign_ = 1 - first_digit;
  }
  BigInteger(bool sign, std::vector<int64_t> digits) : sign_(sign), digits_(digits) {}
  BigInteger(int degree, int fictious) : sign_(true) {  // constructs 10^degree
    for (int i = 0; i < degree / log_base; ++i) {
      digits_.push_back(0);
    }
    std::string remain = "1";
    for (int i = 0; i < degree % log_base; ++i) {
      remain.push_back('0');
    }
    if (fictious != -1) {
      digits_.push_back(std::stoi(remain));
    }
  }
  int64_t& operator[](size_t i) { return digits_[i]; }
  BigInteger operator-() const { return BigInteger(!sign_, digits_); }
  BigInteger& operator+=(const BigInteger& num);
  BigInteger& operator-=(const BigInteger& num);
  BigInteger& operator*=(const BigInteger& num);
  BigInteger& operator/=(const BigInteger& num);
  BigInteger& operator%=(const BigInteger& num);
  BigInteger& operator++();
  BigInteger& operator--();
  BigInteger operator++(int);
  BigInteger operator--(int);
  const int64_t& operator[](size_t i) const { return digits_[i]; }
  std::string toString() const;
  bool sign() const { return sign_; }
  int length() const { return digits_.size(); }
  bool empty() const { return digits_.empty(); }
  BigInteger& clear_zeros();
  BigInteger abs() const { return BigInteger(true, digits_); }
  void add_digit(int64_t digit) { digits_.push_back(digit); }
};

std::string to_base(int digit) {
  std::string new_digit = std::to_string(digit);
  while (new_digit.length() != BigInteger::log_base) {
    new_digit = "0" + new_digit;
  }
  return new_digit;
}

std::string BigInteger::toString() const {
  std::string num = "";
  if (digits_.empty()) {
      num.push_back('0');
      return num;
    }
    if (!sign_) {
      num.push_back('-');
    }
    num += std::to_string(digits_[length() - 1]);
    for (int i = length() - 2; i >= 0; --i) {
      num += to_base(digits_[i]);
    }
    return num;
}
BigInteger& BigInteger::clear_zeros() {
  while (length() > 0 && digits_.back() == 0) {
    digits_.pop_back();
  }
  return *this;
}
BigInteger operator "" _bi(unsigned long long num) { return BigInteger(num); }

bool operator==(const BigInteger& num1, const BigInteger& num2) {
  if (num1.empty() && num2.empty()) {
    return true;
  }
  if (num1.sign() != num2.sign()) {
    return false;
  }
  if (num1.length() != num2.length()) {
    return false;
  }
  for (int i = 0; i < num1.length(); ++i) {
    if (num1[i] != num2[i]) {
      return false;
    }
  }
  return true;
}
bool operator!=(const BigInteger& num1, const BigInteger& num2) {
  return !(num1 == num2);
}
bool operator<(const BigInteger& num1, const BigInteger& num2) {
  if (num1.sign() != num2.sign()) {
    return (!num1.sign() && num2.sign());
  }
  if (num1.length() != num2.length()) {
    return ((num1.length() < num2.length()) == num1.sign());
  }
  for (int i = num1.length() - 1; i >= 0; --i) {
    if (num1[i] == num2[i]) {
      continue;
    }
    return ((num1[i] < num2[i]) == num1.sign());
  }
  return false;
}
bool operator>(const BigInteger& num1, const BigInteger& num2) {
  return num2 < num1;
}
bool operator<=(const BigInteger& num1, const BigInteger& num2) {
  return !(num2 < num1);
}
bool operator>=(const BigInteger& num1, const BigInteger& num2) {
  return !(num1 < num2);
}
bool module_comparison(const BigInteger& num1, const BigInteger& num2) {
  if (num1.length() != num2.length()) {
    return (num1.length() < num2.length());
  }
  for (int i = num1.length() - 1; i >= 0; --i) {
    if (num1[i] == num2[i]) {
      continue;
    }
    return (num1[i] < num2[i]);
  }
  return false;
}

BigInteger& addition(BigInteger& num1, const BigInteger& num2) {  // both num1 and num2 considered to be non-negative
  int64_t carry = 0;
  for (int i = 0; i < std::max(num1.length(), num2.length()); ++i) {
    int64_t digit1 = i < num1.length() ? num1[i] : 0;
    int64_t digit2 = i < num2.length() ? num2[i] : 0;
    if (i < num1.length()) {
      num1[i] = (digit1 + digit2 + carry) % BigInteger::base;
    } else {
      num1.add_digit((digit1 + digit2 + carry) % BigInteger::base);
    }
    carry = (digit1 + digit2 + carry) / BigInteger::base;
  }
  if (carry) {
    num1.add_digit(carry);
  }
  return num1.clear_zeros();
}
BigInteger& subtraction(BigInteger& num1, const BigInteger& num2) {  // --//-- and num1 >= num2
  int64_t carry = 0;
  for (int i = 0; i < std::max(num1.length(), num2.length()); ++i) {
    int64_t digit1 = i < num1.length() ? num1[i] : 0;
    int64_t digit2 = i < num2.length() ? num2[i] : 0;
    int64_t current_digit = digit1 - digit2 - carry;
    int64_t add_base = 0;
    if (current_digit < 0) {
      add_base = BigInteger::base;
      carry = 1;
    } else {
      carry = 0;
    }
    if (i < num1.length()) {
      num1[i] = (current_digit + add_base);
    } else {
      num1.add_digit(current_digit + add_base);
    }
  }
  return num1.clear_zeros();
}
BigInteger operator+(BigInteger num1, const BigInteger& num2) {
  return num1 += num2;
}
BigInteger operator-(BigInteger num1, const BigInteger& num2) {
  return num1 -= num2;
}
BigInteger& BigInteger::operator+=(const BigInteger& num) {
  if (sign() == num.sign()) {
    return  addition(*this, num);;
  }
  if (module_comparison(*this, num)) {
    BigInteger temp = num;
    *this = subtraction(temp, *this);
    return *this;
  }
  return subtraction(*this, num);
}
BigInteger& BigInteger::operator-=(const BigInteger& num) {
  *this += (-num);
  return *this;
}
void BigInteger::increase() {
  int i = 0;
  while (i < length() && digits_[i] == BigInteger::base - 1) {
    std::cout<<1;
    digits_[i] = 0;
    ++i;
  }
  if (i == length()) {
    digits_.push_back(1);
  } else {
    ++digits_[i];
  }
}
void BigInteger::decrease() {
  int i = 0;
  while (i < length() && digits_[i] == 0) {
    digits_[i] = BigInteger::base - 1;
    ++i;
  }
  --digits_[i];
  clear_zeros(); 
}
BigInteger& BigInteger::operator++() {
  if (sign_) {
    increase();
  } else {
    decrease();
  }
  return *this;
}
BigInteger& BigInteger::operator--() {
  if (!sign_) {
    increase();
  } else {
    decrease();
  }
  return *this;
}
BigInteger BigInteger::operator++(int) {
  BigInteger temp(*this);
  ++(*this);
  return temp;
}
BigInteger BigInteger::operator--(int) {
  BigInteger temp(*this);
  --(*this);
  return temp;
}

BigInteger& multiplication(BigInteger& num1, const BigInteger& num2) {  // both num1 and num2 considered to be non-negative
  std::vector<int64_t> product_digits;
  int64_t carry = 0;
  for (int i = 0; i < num1.length() + num2.length() - 1; ++i) {
    int64_t current_digit = 0;
    for (int j = std::min(i, num1.length() - 1); j >= std::max(0, i - num2.length() + 1); --j) {
      current_digit += num1[j] * num2[i - j];
    }
    product_digits.push_back((current_digit + carry) % BigInteger::base);
    carry = (current_digit + carry) / BigInteger::base;
  }
  if (carry) {
    product_digits.push_back(carry);
  }
  num1 = BigInteger(true, product_digits).clear_zeros();
  return num1;
}
BigInteger operator*(BigInteger num1, const BigInteger& num2) {
  num1 *= num2;
  return num1;
}
BigInteger& BigInteger::operator*=(const BigInteger& num) {
  *this = (sign() == num.sign()) ? multiplication(*this, num)
                                : -multiplication(*this, num);
  return *this;
}

std::pair<BigInteger, BigInteger> division(const BigInteger& num1, const BigInteger& num2) {
  std::vector<int64_t> ratio_digits;
  BigInteger temp = 0;
  for (int i = num1.length() - 1; i >= 0; --i) {
    temp *= BigInteger::base;
    temp += num1[i];
    int64_t left = 0, right = BigInteger::base;
    if (temp < num2) {
      ratio_digits.push_back(0);
      continue;
    }
    while (right - left > 1) {
      int64_t m = (left + right) / 2;
      if (m * num2 <= temp) {
        left = m;
      } else {
        right = m;
      }
    }
    temp -= left * num2;
    ratio_digits.push_back(left);
  }
  if (ratio_digits.empty()) {
    return {BigInteger(0), BigInteger(0)};
  }
  std::reverse(ratio_digits.begin(), ratio_digits.end());
  while (ratio_digits.back() == 0) {
    ratio_digits.pop_back();
  }
  return {BigInteger(1, ratio_digits), temp};
}

BigInteger operator/(const BigInteger& num1, const BigInteger& num2) {
  return (num1.sign() == num2.sign()) ? division(num1.abs(), num2.abs()).first
                                      : -division(num1.abs(), num2.abs()).first;
}
BigInteger& BigInteger::operator/=(const BigInteger& num) {
  *this = (*this) / num;
  return *this;
}
BigInteger operator%(const BigInteger& num1, const BigInteger& num2) {
  return num1 - (num1 / num2) * num2;
}
BigInteger& BigInteger::operator%=(const BigInteger& num) {
  *this = *this % num;
  return *this;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& num) {
  if (num.empty()) {
    out << 0;
    return out;
  }
  if (!num.sign()) {
    out << '-';
  }
  out << num[num.length() - 1];
  for (int i = num.length() - 2; i >= 0; --i) {
    out << to_base(num[i]);
  }
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& num) {
  std::string new_num;
  do {
    char digit;
    in >> digit;
    new_num.push_back(digit);
  } while (in.peek() != EOF && !std::isspace(in.peek()));
  num = BigInteger(new_num);
  return in;
}

class Rational {
 private:
  BigInteger numerator_;
  BigInteger denominator_;
  void make_coprime();
 public:
  Rational() : numerator_(0), denominator_(1) {}
  Rational(int64_t num) : numerator_(BigInteger(num)), denominator_(1) {}
  Rational(const BigInteger& num) : numerator_(num), denominator_(1) {}
  Rational(const BigInteger& num, const BigInteger& denom) : numerator_(num), denominator_(denom) {}
  explicit operator double() const { return std::stod((*this).asDecimal(20)); }
  Rational& operator+=(const Rational& num);
  Rational& operator-=(const Rational& num);
  Rational& operator*=(const Rational& num);
  Rational& operator/=(const Rational& num);
  Rational operator-() const { return Rational(-numerator_, denominator_); }
  const BigInteger& num() const {return numerator_; }
  const BigInteger& denom() const {return denominator_; }
  std::string toString() const;
  std::string asDecimal(size_t precision) const;
};


std::string Rational::asDecimal(size_t precision) const {
  BigInteger temp = numerator_ * BigInteger(precision, 0);
  temp /= denominator_;
  std::string decimal = "";
  if (!temp.sign()) {
    decimal.push_back('-');
  }
  if (temp.empty()) {
    return "0";
  }
  std::string digits = (temp.abs()).toString();
  std::reverse(digits.begin(), digits.end());
  while (digits.length() <= precision) {
    digits.push_back('0');
  }
  for (int j = digits.length() - 1; j >= 0; --j) {
    decimal.push_back(digits[j]);
    if (j != 0 && j == static_cast<int>(precision)) {
      decimal.push_back('.');
    }
  }
  return decimal; 
}

std::string Rational::toString() const {
  Rational temp = *this;
  temp.make_coprime();
  std::string num = "";
  if(!temp.numerator_.empty() && temp.numerator_.sign() != temp.denominator_.sign()) {
    num.push_back('-');
  }
  num += (temp.numerator_.abs()).toString();
  if (temp.denominator_ != 1 && temp.denominator_ != -1) {
    num.push_back('/');
    num += (temp.denominator_.abs()).toString();
  }
  return num;
}

std::ostream& operator<<(std::ostream& out, const Rational& num) {
  out << num.asDecimal(5);
  return out;
}
void Rational::make_coprime() {
  BigInteger copy1 = numerator_.abs();
  BigInteger copy2 = denominator_.abs();
  if (copy1 == 1 || copy2 == 1) {
    return;
  }
  while (!copy1.empty() && !copy2.empty()) {
    if (copy1 >= copy2) {
      copy1 -= copy2 * (copy1 / copy2);
    } else {
      copy2 -= copy1 * (copy2 / copy1);
    }
  }
  BigInteger gcd = (!copy1.empty()) ? copy1 : copy2;
  numerator_ /= gcd;
  denominator_ /= gcd;
}
Rational& Rational::operator+=(const Rational& num) {
  numerator_ = numerator_ * num.denominator_
             + denominator_ * num.numerator_;
  denominator_ *= num.denominator_;
  return *this;
}
Rational& Rational::operator-=(const Rational& num) {
  Rational temp = -num;
  *this += temp;
  make_coprime();
  return *this;
}
Rational& Rational::operator*=(const Rational& num) {
  numerator_ *= num.numerator_;
  denominator_ *= num.denominator_;
  return *this;
}
Rational& Rational::operator/=(const Rational& num) {
  numerator_ *= num.denominator_;
  denominator_ *= num.numerator_;
  make_coprime();
  return *this;
}

Rational operator+(Rational num1, const Rational& num2) {
  num1 += num2;
  return num1;
}
Rational operator-(Rational num1, const Rational& num2) {
  num1 -= num2;
  return num1;
}
Rational operator*(Rational num1, const Rational& num2) {
  num1 *= num2;
  return num1;
}
Rational operator/(Rational num1, const Rational& num2) {
  num1 /= num2;
  return num1;
}

bool operator==(const Rational& num1, const Rational& num2) {
  return num1.num() * num2.denom() == num1.denom() * num2.num();
}
bool operator!=(const Rational& num1, const Rational& num2) {
  return !(num1 == num2);
}
bool operator<(const Rational& num1, const Rational& num2) {
  bool sign = (num1.num() * num2.denom() - num1.denom() * num2.num()).sign();
  if (sign) {
    return num1.denom().sign() != num2.denom().sign();
  }
  return num1.denom().sign() == num2.denom().sign();
}
bool operator>(const Rational& num1, const Rational& num2) {
  return num2 < num1;
}
bool operator<=(const Rational& num1, const Rational& num2) {
  return !(num2 < num1);
}
bool operator>=(const Rational& num1, const Rational& num2) {
  return !(num1 < num2);
}
std::istream& operator>>(std::istream& in, Rational& r) {
  std::string s;
  in >> s;
  r = Rational(BigInteger(s));
  return in;
}


template <size_t N, size_t l = 1, size_t r = N>
struct root {
  static const size_t m = (l + r) / 2;
  static const bool bigger = (m * m >= N);
  static const size_t value = root<N, (bigger ? l : m + 1), (bigger ? m : r)>::value;
};

template <size_t N, size_t val>
struct root<N, val, val> {
  static const size_t value = val;
};

template <size_t N>
const static size_t root_v = root<N>::value;


template <size_t N, size_t k>
struct is_prime {
  const static bool value = (N % k != 0) && is_prime<N, k - 1>::value;
};

template <size_t N>
struct is_prime<N, 1> {
  const static bool value = true;
};

template <size_t N>
struct is_prime<N, 0> {
  const static bool value = true;
};

template <size_t N>
const static bool is_prime_v = is_prime<N, root_v<N>>::value;

template <size_t N>
class Residue {
 private:
  size_t remainder;

 public:
  Residue() : remainder(0) {}
  Residue(int n) : remainder((N + (n % static_cast<int>(N))) % static_cast<int>(N)) {}
  explicit operator int() { return remainder; }
  size_t r() const { return remainder; }
  Residue& operator+=(Residue<N> other);
  Residue& operator-=(Residue<N> other);
  Residue& operator*=(Residue<N> other);
  Residue& operator/=(Residue<N> other);
};

template <size_t N>
bool operator==(Residue<N> r1, Residue<N> r2) { return r1.r() == r2.r(); }


template <size_t N>
Residue<N> operator+(Residue<N> r1, Residue<N> r2) { return Residue<N>((r1.r() + r2.r()) % N); }

template <size_t N>
Residue<N> operator-(Residue<N> r1, Residue<N> r2) {
  int diff = static_cast<int>(r1.r()) - static_cast<int>(r2.r());
  return Residue<N>(diff);
}

template <size_t N>
Residue<N> operator*(Residue<N> r1, Residue<N> r2) { return Residue<N>((r1.r() * r2.r()) % N); }

template <size_t N>
Residue<N> operator/(Residue<N> r1, Residue<N> r2) {
  static_assert(is_prime_v<N>, "NOT A FIELD");
  for (size_t i = 0; i < N; ++i) {
    if ((i * r2.r()) % N == r1.r()) {
      return Residue<N>(i);
    }
  }
  return Residue<N>(0);
}

template <size_t N>
Residue<N>& Residue<N>::operator+=(Residue<N> other) {
  *this = *this + other;
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator-=(Residue<N> other) {
  *this = *this - other;
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator*=(Residue<N> other) {
  *this = (*this) * other;
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator/=(Residue<N> other) {
  *this = (*this) / other;
  return *this;
}

template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& r) {
  std::string s;
  in >> s;
  r = Residue<N>(std::stoi(s));
  return in;
}


template <size_t M, size_t N, typename Field=Rational>
class Matrix {
 private:
  std::vector<std::vector<Field>> matrix;
  void swap_columns(size_t i, size_t j);
  void swap_lines(size_t i, size_t j) { std::swap(matrix[i], matrix[j]); }
  std::pair<size_t, size_t> find_non_zero(size_t i1, size_t j1, size_t i2, size_t j2) const;
  std::pair<Matrix<M, N, Field>, int> make_Gauss() const;

 public:
  static const size_t INF = 1e9;
  Matrix() { matrix.resize(M, std::vector<Field>(N)); }
  Matrix(std::initializer_list<std::initializer_list<Field>> mtx) : Matrix() {
    int i = 0;
    for (std::initializer_list<Field> line : mtx) {
      int j = 0;
      for (Field element : line) {
        matrix[i][j] = element;
        ++j;
      }
      ++i;
    }
  }
  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other);
  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other);
  Matrix<M, N, Field>& operator*=(const Matrix<M, N, Field>& other);

  const std::vector<Field>& operator[](size_t row) const { return matrix[row]; }
  std::vector<Field>& operator[](size_t row) { return matrix[row]; }
  Field trace() const;
  Field det() const;
  size_t rank() const;
  Matrix<N, M, Field> transposed() const;
  void invert();
  Matrix<M, N, Field> inverted() const;
  std::array<Field, N> getRow(size_t i) const;
  std::array<Field, M> getColumn(size_t i) const;
};

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Field f, const Matrix<M, N, Field>& matrix) {
  Matrix<M, N, Field> temp = matrix;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      temp[i][j] *= f;
    }
  }
  return temp;
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::swap_columns(size_t i, size_t j) {
  for (int k = 0; k < M; ++k) {
    std::swap(matrix[k][i], matrix[k][j]);
  }
}

template <size_t M, size_t N, typename Field>
bool operator==(const Matrix<M, N, Field>& m1, const Matrix<M, N, Field>& m2) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (!(m1[i][j] == m2[i][j])) {
        return false;
      }
    }
  }
  return true;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Matrix<M, N, Field>& other) {
  static_assert(M == N);
  Matrix<M, N, Field> result = *this;
  *this = *this * other;
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field>& other) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      (*this)[i][j] += other[i][j];
    }
  }
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field>& other) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      (*this)[i][j] -= other[i][j];
    }
  }
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(Matrix<M, N, Field> m1, const Matrix<M, N, Field>& m2) {
  return m1 += m2; 
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(Matrix<M, N, Field> m1, const Matrix<M, N, Field>& m2) {
  return m1 -= m2; 
}

template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& m1, const Matrix<N, K, Field>& m2) {
  Matrix<M, K, Field> product_matrix;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      Field sum(0);
      for (size_t k = 0; k < N; ++k) {
        sum += m1[i][k] * m2[k][j]; 
      }
      product_matrix[i][j] = sum;
    }
  }
  return product_matrix;
}

template <size_t M, size_t N, typename Field>
std::pair<Matrix<M, N, Field>, int> Matrix<M, N, Field>::make_Gauss() const {  // returns matrix after Gauss transformation and
  Matrix<M, N, Field> new_matrix = *this;                                      //  coefficient +-1(parity of lines transpositions)
  size_t i = 0, j = 0;
  int coeff = 1;
  while (i < M && j < N) {
    auto non_zero = new_matrix.find_non_zero(i, j, M, N);
    if (non_zero == std::make_pair(INF, INF)) {
      return {new_matrix, coeff};
    }
    new_matrix.swap_lines(i, non_zero.first);
    if (i != non_zero.first) {
      coeff *= -1;
    }
    for (size_t down_i = i + 1; down_i < M; ++down_i) {
      Field to_multiply = new_matrix[down_i][non_zero.second] / new_matrix[i][non_zero.second];
      for (size_t down_j = 0; down_j < N; ++down_j) {
        new_matrix[down_i][down_j] -= to_multiply * new_matrix[i][down_j];
      }
    }
    ++i;
    ++j;
  }
  return {new_matrix, coeff};
}

template <size_t M, size_t N, typename Field>
std::pair<size_t, size_t> Matrix<M, N, Field>::find_non_zero(size_t i1, size_t j1, size_t i2, size_t j2) const {  // returns the most left (then the most upper) non-zero element
  for (size_t j = j1; j < j2; ++j) {
    for (size_t i = i1; i < i2; ++i) {
      if (!(matrix[i][j] == Field(0))) {
        return {i, j};
      }
    }
  }
  return {INF, INF};
}

template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
  static_assert(M == N);
  Field trc = matrix[0][0];
  for (size_t i = 1; i < M; ++i) {
    trc += matrix[i][i];
  }
  return trc;
}

template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
  static_assert(M == N);
  auto ladder_matrix = make_Gauss();
  Field mult = (ladder_matrix.first)[0][0];
  for (size_t i = 1; i < M; ++i) {
    mult *= (ladder_matrix.first)[i][i];
  }
  return mult * static_cast<Field>(ladder_matrix.second);
}

template<size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
  auto ladder_matrix = make_Gauss().first;
  size_t rnk = 0;
  for (size_t i = 0; i < std::min(M, N); ++i) {
    bool non_zeros = false;
    for (size_t j = 0; j < N; ++j) {
      if (!(ladder_matrix[i][j] == Field(0))) {
        non_zeros = true;
      }
    }
    if (non_zeros) {
      rnk++;
    } else {
      return rnk;
    }
  }
  return rnk;
}

template<size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
  Matrix<N, M, Field> transposed_matrix;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      transposed_matrix[j][i] = matrix[i][j];
    }
  }
  return transposed_matrix;
}

template<size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
  static_assert(M == N);
  Matrix<M, N, Field> new_matrix;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      new_matrix[i][j] = (i == j ? Field(1) : Field(0));
    }
  }
  size_t i = 0, j = 0;
  while (i < M && j < N) {
    auto non_zero = find_non_zero(i, j, M, N);
    swap_lines(i, non_zero.first);
    new_matrix.swap_lines(i, non_zero.first);
    Field to_divide = matrix[i][non_zero.second];
    for (size_t cur_j = 0; cur_j < N; ++cur_j) {
      matrix[i][cur_j] /= to_divide;
      new_matrix[i][cur_j] /= to_divide;
    }
    for (size_t down_i = 0; down_i < M; ++down_i) {
      if (down_i == i) {
        continue;
      }
      Field to_multiply = matrix[down_i][non_zero.second] / matrix[i][non_zero.second];
      for (size_t down_j = 0; down_j < N; ++down_j) {
        matrix[down_i][down_j] -= to_multiply * matrix[i][down_j];
        new_matrix[down_i][down_j] -= to_multiply * new_matrix[i][down_j];
      }
    }
    ++i;
    ++j;
  }
  *this = new_matrix;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
  Matrix<M, N, Field> result = *this;
  result.invert();
  return result;
}

template<size_t M, size_t N, typename Field>
std::array<Field, M> Matrix<M, N, Field>::getColumn(size_t i) const {
  std::array<Field, M> column;
  for (size_t j = 0; j < M; ++j) {
    column[j] = matrix[j][i];
  }
  return column;
}

template<size_t M, size_t N, typename Field>
std::array<Field, N> Matrix<M, N, Field>::getRow(size_t i) const {
  std::array<Field, N> row;
  for (size_t j = 0; j < N; ++j) {
    row[j] = matrix[i][j];
  }
  return row;
}