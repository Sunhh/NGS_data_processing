package mathSunhh; 
# A package of math/utility sub-routines. The implementation is split across
# mathSunhh/*.pm parts (each declares 'package mathSunhh;') for navigability;
# this file is the loader. All subs live in the single 'mathSunhh' namespace.
# It is function-only: every sub is called as mathSunhh::foo(...); no constructor/objects.
use strict; 
use warnings; 
use Statistics::Descriptive; 
use Scalar::Util qw(looks_like_number blessed);
use Exporter qw(import);
our @EXPORT = qw(ins_calc);
our @EXPORT_OK = qw(parseCol p_adjust_BH chisqrprob);
use LogInforSunhh; 

############################################################
#  Load the split parts (each declares 'package mathSunhh;').
############################################################
require mathSunhh::Param; 
require mathSunhh::Stats; 
require mathSunhh::ArrayComb; 
require mathSunhh::ObjNum; 
require mathSunhh::Windows; 
require mathSunhh::Interval; 
require mathSunhh::Group; 
require mathSunhh::Encode; 
require mathSunhh::Color; 
require mathSunhh::IndexDB; 


1; 
