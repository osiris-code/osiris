# Introduction

Osiris simulations are set up using an input file (also referred to as
input deck throughout this document) with the format specified here.
When Osiris is launched it will look for the specified input file (or a
file named *os-stdin* in its run directory if none was supplied) and
read the simulation parameters from it. If this file is not present the
code stop with an error message. Note that when running on most
platforms (like the eXpp cluster, for example) the startup scripts will
take care of copying and renaming your input file. However, if you are
launching osiris manually, you must do this yourself.

## Command Line Options

The Osiris binaries will accept the following command line options:

```text
Usage: osiris [-t] [-r] [-w work_dir] [input_file]
 
Arguments:
 input_file    -   Input file to use. If not specified defaults to "os-stdin"

Options:

 -t            -   Test only. If set osiris will only test the validity of the specified file, and all other options are ignored.
 -r            -   Restart. If set it will force osiris to restart from checkpoint data.
 -w work_dir   -   Work directory. If set osiris will change to this directory before starting the simulation.
```

Alternatively you can use the following environment variables to control
Osiris behavior:

* `OSIRIS_TEST` - Setting this to any value is the same as setting the `-t` option above.
* `OSIRIS_RESTART` - Setting this to any value is the same as setting the `-r` option above.
* `OSIRIS_WDIR` - Set this to work_dir as in the `-w` option above.

The environment variables override the command line options.

## Input File Format

The Osiris input files are a sequence of sections specifying the several
parameters for the simulation to be performed. The format for these
sections is the following:

```text
section-name
{
 data
} 
```

The section data is specified using the normal Fortran namelist
conventions i.e. "var = value," (don't forget the trailing comma). Here
are two examples of sections on osiris input decks:

```text
! grid parameters 
grid 
{
    nx_p(1:2) = 4096, 512,      ! use a 4096 x 512 cell grid
    coordinates = "cartesian",  ! with cartesian coordinates
}
   
! empty section, use default parameters
current 
{}
```

When any variable of the data section is not specified it default to a
pre-defined value. To find out what these values check the Osiris Input
section below. The opening and closing brackets can be placed in the
same line as the section-name.

### Comments

Any text following an exclamation mark "!" is ignored until the end of
line is reached. Blank lines are also ignored. Examples of comments can
be seen in the above section

## Osiris Input

An Osiris input file is made out of the following sections, placed in
the order defined below. Follow the links for a detailed description of
each section. Sections marked with (optional) are not required for an
osiris run and can be left out of the input file. In this situation the
code will issue a warning message stating that it didn't find the
section and is using the defaults. The order of the sections CAN NOT be
changed and doing so will result in a run-time error.

### General Simulation Parameters

* [simulation](Simulation.md) - (optional) -  global simulation parameters
* [node_conf](Node_Configuration.md) - node configuration, and periodic boundary settings
* [grid](Grid.md) - grid and coordinate system
* [time_step](Time_Step.md) - time step and dump frequency
* [restart](Restart.md) - (optional) - restart settings
* [space](Space.md) - spacial dimensions and moving window
* [time](Time.md) - time limits

### Electro-Magnetic Fields

* [el_mag_fld](Electro-Magnetic_Field.md) - (optional) - solver type, initial/external fields, subcycling and
  smoothing
* [emf_bound](Electro-Magnetic_Field_Boundaries.md) - electro-magnetic field boundaries
* [pgc](Ponderomotive_Guiding_Center_Laser.md) - (optional) - ponderomotive guiding center laser
* [smooth](Smooth.md) - (optional) - electro-magnetic field smoothing
* [diag_emf](Electro-Magnetic_Field_Diagnostics.md) - (optional) - electo-magnetic field diagnostics

### Particles

* [particles](Particles.md) - number of species, cathodes, neutrals and neutral_mov_ions in the simulation and
  global particle diagnostics.

The following items are repeated for every species

* [species](Species.md) - information on the
  particle species
* [udist](Momentum_Distribution.md) -
  (optional) - momentum distribution for the species
* [profile](Profile.md) - (optional) - density
  profile information for the species
* [spe_bound](Species_Boundary.md) - boundary
  conditions for the species
* [diag_species](Species_Diagnostics.md) -
  (optional) - species diagnostics

The following items are repeated for every cathode

* [cathode](Cathode.md) - information on the
  cathode
* [species](Species.md) - information on the
  particle species to inject electrons from the cathode
* [udist](Momentum_Distribution.md) -
  (optional) - momentum distribution for the species
* [spe_bound](Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral

* [neutral](Neutrals.md) - information on
  ionization gas and ion model
* [profile](Profile.md) - (optional) - density
  profile information for the neutral gas
* [diag_neutral](Neutrals_Diagnostics.md) -
  (optional) - neutral gas diagnostics
* [species](Species.md) - information on the
  particle species to inject electrons from ionization
* [spe_bound](Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral_mov_ions

* [neutral_mov_ions](Neutrals_with_Moving_Ions.md) -
  information on ionization gas and ion model
* [profile](Profile.md) - (optional) - density
  profile information for the neutral gas
* [diag_neutral](Neutrals_Diagnostics.md) -
  (optional) - neutral gas diagnostics
* [species](Species.md) - information on the
  particle species to inject electrons from ionization
* [spe_bound](Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons
* [species](Species.md) - information on the
  particle species to inject ions from ionization
* [spe_bound](Species_Boundary.md) - boundary
  conditions for the species holding injected ions
* [diag_species](Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected ions

The following item is specified after all the species, cathode, neutral,
and neutral_mov_ions information.

* [collisions](Collisions.md) - (optional)
  information on the binary collisions between particles

### Laser Pulses

The following item is repeated for every laser pulse to be launched.
There is no section specifying the number of laser pulses; this is
determined automatically from the number of zpulse sections in the input
deck.

* [zpulse](Laser_Pulse.md) - (optional) -
  laser pulse parameters

### Electrical Current

* [current](Electric_Current.md) -
  (optional) - dummy section for separating the electrical current
  pfield section
* [smooth](Smooth.md) - (optional) - current
  smoothing information
* [diag_current](Electric_Current_Diagnostics.md)-
  (optional) - electrical current diagnostics

### Antennas

* [antenna_array](Antenna_Array.md)-
  (optional) - number of antennas to use

The following item is repeated for every antenna

* [antenna](Antenna.md) - information on the
  antenna to use

## Analytical Function Parser

Osiris includes an analytical function parser, so that the user can
specify parameters for the simulation in the form of an analytical
expression that is evaluated at run time.

This function will be compiled into a pseudo-code that evaluates quite
fast. However the user should bear in mind that this will always be
slower to execute than a Fortran function that is compiled into Osiris.
Furthermore, you should note that compilation is done without any kind
of optimization, so that for better performance all unnecessary
calculations should be removed from the expression (i.e. "6.4^2" should
be written as "40.96").

All calculations are done using double-precision floating point
arithmetic, i.e., using the default Fortran type `real(kind(1.0d0))`. This
data will from now on be referred to as numeric. As an extension these
routines a logical data type is also included, implemented in a manner
similar to the C programming language. A value of 0.0 is considered to
be false, and all other values are considered to be true. This data will
from now on be referred to as logical, and is relevant to the `if`
function described below.

The following operators are defined for this function parser, shown in
decreasing order of relative precedence.

Relative precedence of operators (in increasing order)

| Type | Operator |
| :---       | :-: |
| logical    | `&&`, `\|\|` |
| logical    | `==`, `!=`, `>=`, `<=`, `>`, `<` |
| numeric    | monadic `+` or `-` |
| numeric    | dyadic `+` or `-` |
| numeric    | `*` or `/` |
| numeric    | `^` (exponentiation) |

Parenthesis "( )" can (must) be used to change the priorities of
operators within an expression. These operators implement the following
operations:

* `^` - Exponentiation, "a^b" represents a to the power of b. (Note that b
  can be real here but is cast to an int; use the pow() function for
  non-integer exponents.)
* `+`, `-` (monadic) - "+a" returns a, and "-a" returns -a.
* `+`, `-` (dyadic) - Sum, subtraction, "a+b" returns a plus b and "a-b"
  returns a minus b.
* `*`,`/` - Multiplication, division, "a\*b" returns a times b and "a/b"
  returns a divided by b (real division). Integer division is not
  implemented but can be done using the floor function described below.
* `==` - Equal, "a==b" returns 1.0 if a is equal to b and 0.0 otherwise.
* `!=` - Different, "a!=b" returns 1.0 if a is different from b and 0.0
  otherwise.
* `>` - Greater than, "a\>b" returns 1.0 if a is greater than b and 0.0
  otherwise.
* `<` - Smaller than, "a\<b" returns 1.0 if a is smaller than b and 0.0
  otherwise.
* `>=` - Greater than or equal, "a\>=b" returns 1.0 if a is greater than
  b or equal to b and 0.0 otherwise.
* `<=` - Smaller than or equal, "a\<=b" returns 1.0 if a is smaller than
  b or equal to b and 0.0 otherwise.
* `&&` - Logical intersection, "a && b" returns 1.0 if a and b are true
  and 0.0 otherwise. a and b are treated as logical values.
* `||` - Logical union, "a \|\| b" returns 1.0 if a or b are true and
  0.0 otherwise. a and b are treated as logical values.

Here are some examples of simple expressions:

"x1 + (2.0*x2 - 1.0)^2"
"-x3/2.0 + x1^2"

This parser also allows for the use of an "if" function for the
conditional evaluation of expressions. The syntax for this function is
the following:

* `if( test, A, B )`

`test` is evaluated as a logical value. If test is true then `A` is
returned, otherwise `B` is returned. `A` and `B` can be any valid expressions,
including other if functions, thus allowing for nested if structures.
Here are some examples:

```text
"if( x1 > 5.0, x1-5.0, 0.0)"
"if( (x1 > 4.0) && (x1<5.0), 1.0, 0.0)"
```

Finally the following mathematical functions are also implemented. The
general syntax is:

* `func(x)`, `func(x,y)`, `func(x,y,z)` - for single, double or triple parameter functions

Note that for relevant trigonometric functions the parameter/result is
in radians.

* `abs(x)` - Absolute value of x.
* `sin(x)` - Sine of x.
* `cos(x)` - Cosine of x.
* `tan(x)` - Tangent of x.
* `exp(x)` - Exponential function i.e. e^x.
* `log10(x)` - Base 10 logarithm of x.
* `log(x)` - Natural (Base e) logarithm of x.
* `asin(x)` - Arc Sine of x.
* `acos(x)` - Arc Cosine of x.
* `atan2(x,y)` - Arc Tangent of y/x, taking into account which quadrant
  the point (x,y) is in.
* `atan(x)` - Arc Tangent of x.
* `sqrt(x)` - Square root of x.
* `not(x)` - Logical not. x is evaluated as a logical expression and the
  complement is returned.
* `pow(x,y)` - Power, returns x^y.
* `int(x)` - Integer, converts x to integer truncating towards 0.
* `nint(x)` - Nearest integer, converts x to the nearest integer.
* `ceiling(x)` - Ceiling, converts x to the least integer that is \>= x.
* `floor(x)` - Floor, converts x to the greatest integer that is \<= x.
* `modulo(x,y)` - Modulo, returns the remainder of the integer division,
  i.e., `x - floor(x/y)*y`
* `rect(x)` - Rect function, returns 1.0 for 0.5\<= x \<= 0.5 and 0.0
  otherwise.
* `step(x)` - Step function, returns 1.0 for x \>= 0 and 0.0 otherwise
* `min3(x,y,z)` - Minimum function, returns the minimum value between x, y
  and z
* `min(x,y)` - Minimum function, returns the minimum value between x and y
* `max3(x,y,z)` - Maximum function, returns the minimum value between x, y
  and z
* `max(x,y)` - Maximum function, returns the minimum value between x and y

Here's an example defining a pac-man shaped density in 2D:

```text
math_func_expr = "if((x1-6.4)^2+(x2-6.4)^2<6.4^2,
                                  if( (abs(x2-6.4) < 4.5-x1),0.,1.),
                                  0.0)",
```

## Grid Diagnostics

Starting with revision 357 diagnostics of grid quantities in OSIRIS
share a common interface that allows the user to do several different
types of grid diagnostics besides the simple dumping of the grid
quantity to disk. The user may choose to do data reduction operations (
such as spatial averaging and lineouts ) and also to perform time
averaging operations to filter out higher (time) frequency components.
For a full description of the syntax and available diagnostic types see
the [grid diagnostics section](Grid_Diagnostics.md).

## Input File Examples

The following are example input files of full OSIRIS simulations. The
files will list at the top the version in which they where verified.

* [epc001.2d](MISSING) - 2D Electron-Positiron Weibel
  instability simulation. This input file models the collision of an
  electron cloud with a positron cloud, moving out of the plane with
  $u_{fl} = 0.6 \, \rm{c}$, and thermal spread $u_{th} = 0.1 \, \rm{c}$.

<!-- -->

* [pwfa.2dcyl](MISSING) - 2D cylindrical plasma
  wakefield accelerator, using an electron beam driver with
  $\gamma = 3 \times 10^5$. This also illustrates the use of the
  *n_accelerate* parameter. See the documentation on
  the [udist](Momentum_Distribution.md) section
  for details.

<!-- -->

* [lwfa.ion.2d](MISSING) - 2D Laser wakefield
  accelerator, using a neutral hydrogen gas target, and a laser beam
  driver with an a0 of 10. This illustrates the use of the ionization
  and zpulse modules.
