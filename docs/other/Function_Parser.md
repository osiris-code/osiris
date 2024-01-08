---
layout: single
classes: wide
title: Analytical Function Parser
permalink: /other/function_parser
usemathjax: true

sidebar:
  nav: "other"
---

OSIRIS includes an analytical function parser so that the user can specify parameters for the simulation in the form of an analytical expression that is evaluated at run time.

This function will be compiled into a pseudo-code that evaluates quite fast. However, the user should bear in mind that this will always be slower to execute than a Fortran function that is compiled into OSIRIS. Furthermore, you should note that compilation is done without any kind of optimization, so that for better performance all unnecessary
calculations should be removed from the expression (i.e. "6.4^2" should be written as "40.96").

All calculations are done using double-precision floating-point arithmetic, i.e., using the default Fortran type `real(kind(1.0d0))`. This data will from now on be referred to as numeric. As an extension of these routines, a logical data type is also included, implemented like in the C programming language. A value of 0.0 is considered to be false, and all other values are considered to be true. This data will from now on be referred to as logical and is relevant to the `if` function described below.

The following operators are defined for this function parser, shown in decreasing order of relative precedence.

**Relative precedence of operators (in increasing order)**

| Type | Operator |
| :---       | :-: |
| logical    | `&&`, `||` |
| logical    | `==`, `!=`, `>=`, `<=`, `>`, `<` |
| numeric    | monadic `+` or `-` |
| numeric    | dyadic `+` or `-` |
| numeric    | `*` or `/` |
| numeric    | `^` (exponentiation) |

Parenthesis "( )" can (must) be used to change the priorities of operators within an expression. These operators implement the following operations:

* `^` - Exponentiation, `a^b` represents a to the power of b. (Note that b can be real here but is cast to an int; use the pow() function for non-integer exponents.)
* `+`, `-` (monadic) - `+a` returns a, and `-a` returns -a.
* `+`, `-` (dyadic) - Sum, subtraction, `a+b` returns a plus b and `a-b` returns a minus b.
* `*`,`/` - Multiplication, division, `a*b` returns a times b and `a/b` returns a divided by b (real division). Integer division is not implemented but can be done using the floor function described below.
* `==` - Equal, `a==b` returns 1.0 if a is equal to b and 0.0 otherwise.
* `!=` - Different, `a!=b` returns 1.0 if a is different from b and 0.0 otherwise.
* `>` - Greater than, `a>b` returns 1.0 if a is greater than b and 0.0 otherwise.
* `<` - Smaller than, `a<b` returns 1.0 if a is smaller than b and 0.0 otherwise.
* `>=` - Greater than or equal, `a\>=b` returns 1.0 if a is greater than b or equal to b and 0.0 otherwise.
* `<=` - Smaller than or equal, `a\<=b` returns 1.0 if a is smaller than b or equal to b and 0.0 otherwise.
* `&&` - Logical intersection, `a && b` returns 1.0 if a and b are true and 0.0 otherwise. a and b are treated as logical values.
* `||` - Logical union, `a \|\| b` returns 1.0 if a or b is true and 0.0 otherwise. a and b are treated as logical values.

Here are some examples of simple expressions:

```text
"x1 + (2.0*x2 - 1.0)^2"
"-x3/2.0 + x1^2"
```

This parser also allows for the use of an "if" function for the conditional evaluation of expressions. The syntax for this function is the following:

* `if( test, A, B )`

`test` is evaluated as a logical value. If `test` is true then `A` is returned, otherwise, `B` is returned. `A` and `B` can be any valid expressions, including other if functions, thus allowing for nested if structures. Here are some examples:

```text
"if( x1 > 5.0, x1-5.0, 0.0)"
"if( (x1 > 4.0) && (x1<5.0), 1.0, 0.0)"
```

Finally, the following mathematical functions are also implemented. The general syntax is:

* `func(x)`, `func(x,y)`, `func(x,y,z)` - for single, double or triple parameter functions

Note that for relevant trigonometric functions the parameter/result is in radians.

* `abs(x)` - Absolute value of x.
* `sin(x)` - Sine of x.
* `cos(x)` - Cosine of x.
* `tan(x)` - Tangent of x.
* `exp(x)` - Exponential function i.e. e^x.
* `log10(x)` - Base 10 logarithm of x.
* `log(x)` - Natural (Base e) logarithm of x.
* `asin(x)` - Arc Sine of x.
* `acos(x)` - Arc Cosine of x.
* `atan2(y,x)` - Arc Tangent of y/x, taking into account which quadrant the point (x,y) is in.
* `atan(x)` - Arc Tangent of x.
* `sqrt(x)` - Square root of x.
* `not(x)` - Logical not. x is evaluated as a logical expression and the complement is returned.
* `pow(x,y)` - Power, returns x^y.
* `int(x)` - Integer, converts x to integer truncating towards 0.
* `nint(x)` - Nearest integer, converts x to the nearest integer.
* `ceiling(x)` - Ceiling, converts x to the least integer that is \>= x.
* `floor(x)` - Floor, converts x to the greatest integer that is \<= x.
* `modulo(x,y)` - Modulo, returns the remainder of the integer division, i.e., `x - floor(x/y)*y`
* `rect(x)` - Rect function, returns 1.0 for 0.5\<= x \<= 0.5 and 0.0 otherwise.
* `step(x)` - Step function, returns 1.0 for x \>= 0 and 0.0 otherwise `min3(x,y,z)` - Minimum function, returns the minimum value between x, y and z
* `min(x,y)` - Minimum function, returns the minimum value between x and y
* `max3(x,y,z)` - Maximum function, returns the minimum value between x, y
  and z
* `max(x,y)` - Maximum function, returns the minimum value between x and y

Here's an example defining a pac-man shaped density in 2D:

```text
math_func_expr = "if( (x1-6.4)^2+(x2-6.4)^2<6.4^2,
                     if( (abs(x2-6.4) < 4.5-x1), 0., 1. ),
                     0.0 )",
```
