��7
��
B
AssignVariableOp
resource
value"dtype"
dtypetype�
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
�
Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�
?
Mul
x"T
y"T
z"T"
Ttype:
2	�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
<
Selu
features"T
activations"T"
Ttype:
2
H
ShardedFilename
basename	
shard

num_shards
filename
3
Square
x"T
y"T"
Ttype:
2
	
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
�
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
<
Sub
x"T
y"T
z"T"
Ttype:
2	
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.7.02v2.7.0-rc1-69-gc256c071bb28��5
d
VariableVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable
]
Variable/Read/ReadVariableOpReadVariableOpVariable*
_output_shapes
: *
dtype0
h

Variable_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_1
a
Variable_1/Read/ReadVariableOpReadVariableOp
Variable_1*
_output_shapes
: *
dtype0
h

Variable_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_2
a
Variable_2/Read/ReadVariableOpReadVariableOp
Variable_2*
_output_shapes
: *
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
{
dense_16/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�* 
shared_namedense_16/kernel
t
#dense_16/kernel/Read/ReadVariableOpReadVariableOpdense_16/kernel*
_output_shapes
:	�*
dtype0
s
dense_16/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_16/bias
l
!dense_16/bias/Read/ReadVariableOpReadVariableOpdense_16/bias*
_output_shapes	
:�*
dtype0
|
dense_17/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��* 
shared_namedense_17/kernel
u
#dense_17/kernel/Read/ReadVariableOpReadVariableOpdense_17/kernel* 
_output_shapes
:
��*
dtype0
s
dense_17/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_17/bias
l
!dense_17/bias/Read/ReadVariableOpReadVariableOpdense_17/bias*
_output_shapes	
:�*
dtype0
|
dense_18/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��* 
shared_namedense_18/kernel
u
#dense_18/kernel/Read/ReadVariableOpReadVariableOpdense_18/kernel* 
_output_shapes
:
��*
dtype0
s
dense_18/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_18/bias
l
!dense_18/bias/Read/ReadVariableOpReadVariableOpdense_18/bias*
_output_shapes	
:�*
dtype0
|
dense_19/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��* 
shared_namedense_19/kernel
u
#dense_19/kernel/Read/ReadVariableOpReadVariableOpdense_19/kernel* 
_output_shapes
:
��*
dtype0
s
dense_19/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_19/bias
l
!dense_19/bias/Read/ReadVariableOpReadVariableOpdense_19/bias*
_output_shapes	
:�*
dtype0
|
dense_20/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��* 
shared_namedense_20/kernel
u
#dense_20/kernel/Read/ReadVariableOpReadVariableOpdense_20/kernel* 
_output_shapes
:
��*
dtype0
s
dense_20/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_20/bias
l
!dense_20/bias/Read/ReadVariableOpReadVariableOpdense_20/bias*
_output_shapes	
:�*
dtype0
|
dense_21/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��* 
shared_namedense_21/kernel
u
#dense_21/kernel/Read/ReadVariableOpReadVariableOpdense_21/kernel* 
_output_shapes
:
��*
dtype0
s
dense_21/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_21/bias
l
!dense_21/bias/Read/ReadVariableOpReadVariableOpdense_21/bias*
_output_shapes	
:�*
dtype0
{
dense_22/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�* 
shared_namedense_22/kernel
t
#dense_22/kernel/Read/ReadVariableOpReadVariableOpdense_22/kernel*
_output_shapes
:	�*
dtype0
r
dense_22/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_22/bias
k
!dense_22/bias/Read/ReadVariableOpReadVariableOpdense_22/bias*
_output_shapes
:*
dtype0
z
dense_23/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:* 
shared_namedense_23/kernel
s
#dense_23/kernel/Read/ReadVariableOpReadVariableOpdense_23/kernel*
_output_shapes

:*
dtype0
r
dense_23/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_23/bias
k
!dense_23/bias/Read/ReadVariableOpReadVariableOpdense_23/bias*
_output_shapes
:*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
�
Adam/dense_16/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*'
shared_nameAdam/dense_16/kernel/m
�
*Adam/dense_16/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_16/kernel/m*
_output_shapes
:	�*
dtype0
�
Adam/dense_16/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_16/bias/m
z
(Adam/dense_16/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_16/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_17/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_17/kernel/m
�
*Adam/dense_17/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_17/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/dense_17/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_17/bias/m
z
(Adam/dense_17/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_17/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_18/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_18/kernel/m
�
*Adam/dense_18/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_18/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/dense_18/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_18/bias/m
z
(Adam/dense_18/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_18/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_19/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_19/kernel/m
�
*Adam/dense_19/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_19/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/dense_19/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_19/bias/m
z
(Adam/dense_19/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_19/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_20/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_20/kernel/m
�
*Adam/dense_20/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_20/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/dense_20/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_20/bias/m
z
(Adam/dense_20/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_20/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_21/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_21/kernel/m
�
*Adam/dense_21/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_21/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/dense_21/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_21/bias/m
z
(Adam/dense_21/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_21/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_22/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*'
shared_nameAdam/dense_22/kernel/m
�
*Adam/dense_22/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_22/kernel/m*
_output_shapes
:	�*
dtype0
�
Adam/dense_22/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/dense_22/bias/m
y
(Adam/dense_22/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_22/bias/m*
_output_shapes
:*
dtype0
�
Adam/dense_23/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*'
shared_nameAdam/dense_23/kernel/m
�
*Adam/dense_23/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_23/kernel/m*
_output_shapes

:*
dtype0
�
Adam/dense_23/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/dense_23/bias/m
y
(Adam/dense_23/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_23/bias/m*
_output_shapes
:*
dtype0
�
Adam/dense_16/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*'
shared_nameAdam/dense_16/kernel/v
�
*Adam/dense_16/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_16/kernel/v*
_output_shapes
:	�*
dtype0
�
Adam/dense_16/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_16/bias/v
z
(Adam/dense_16/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_16/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_17/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_17/kernel/v
�
*Adam/dense_17/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_17/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/dense_17/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_17/bias/v
z
(Adam/dense_17/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_17/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_18/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_18/kernel/v
�
*Adam/dense_18/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_18/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/dense_18/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_18/bias/v
z
(Adam/dense_18/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_18/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_19/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_19/kernel/v
�
*Adam/dense_19/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_19/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/dense_19/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_19/bias/v
z
(Adam/dense_19/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_19/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_20/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_20/kernel/v
�
*Adam/dense_20/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_20/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/dense_20/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_20/bias/v
z
(Adam/dense_20/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_20/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_21/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*'
shared_nameAdam/dense_21/kernel/v
�
*Adam/dense_21/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_21/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/dense_21/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*%
shared_nameAdam/dense_21/bias/v
z
(Adam/dense_21/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_21/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_22/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*'
shared_nameAdam/dense_22/kernel/v
�
*Adam/dense_22/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_22/kernel/v*
_output_shapes
:	�*
dtype0
�
Adam/dense_22/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/dense_22/bias/v
y
(Adam/dense_22/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_22/bias/v*
_output_shapes
:*
dtype0
�
Adam/dense_23/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*'
shared_nameAdam/dense_23/kernel/v
�
*Adam/dense_23/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_23/kernel/v*
_output_shapes

:*
dtype0
�
Adam/dense_23/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/dense_23/bias/v
y
(Adam/dense_23/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_23/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
�R
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�Q
value�QB�Q B�Q
�
a1
a2
a3
encoder
decoder
	optimizer
	variables
trainable_variables
	regularization_losses

	keras_api

signatures
;9
VARIABLE_VALUEVariablea1/.ATTRIBUTES/VARIABLE_VALUE
=;
VARIABLE_VALUE
Variable_1a2/.ATTRIBUTES/VARIABLE_VALUE
=;
VARIABLE_VALUE
Variable_2a3/.ATTRIBUTES/VARIABLE_VALUE
�
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
layer_with_weights-4
layer-4
layer_with_weights-5
layer-5
layer_with_weights-6
layer-6
	variables
trainable_variables
regularization_losses
	keras_api
y
layer_with_weights-0
layer-0
	variables
trainable_variables
regularization_losses
	keras_api
�
iter

beta_1

beta_2
	decay
 learning_rate!m�"m�#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�-m�.m�/m�0m�!v�"v�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�-v�.v�/v�0v�
�
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015
16
17
18
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015
 
�
1non_trainable_variables

2layers
3metrics
4layer_regularization_losses
5layer_metrics
	variables
trainable_variables
	regularization_losses
 
h

!kernel
"bias
6	variables
7trainable_variables
8regularization_losses
9	keras_api
h

#kernel
$bias
:	variables
;trainable_variables
<regularization_losses
=	keras_api
h

%kernel
&bias
>	variables
?trainable_variables
@regularization_losses
A	keras_api
h

'kernel
(bias
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
h

)kernel
*bias
F	variables
Gtrainable_variables
Hregularization_losses
I	keras_api
h

+kernel
,bias
J	variables
Ktrainable_variables
Lregularization_losses
M	keras_api
h

-kernel
.bias
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
f
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
f
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
 
�
Rnon_trainable_variables

Slayers
Tmetrics
Ulayer_regularization_losses
Vlayer_metrics
	variables
trainable_variables
regularization_losses
h

/kernel
0bias
W	variables
Xtrainable_variables
Yregularization_losses
Z	keras_api

/0
01

/0
01
 
�
[non_trainable_variables

\layers
]metrics
^layer_regularization_losses
_layer_metrics
	variables
trainable_variables
regularization_losses
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEdense_16/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEdense_16/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEdense_17/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEdense_17/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEdense_18/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEdense_18/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEdense_19/kernel&variables/6/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEdense_19/bias&variables/7/.ATTRIBUTES/VARIABLE_VALUE
KI
VARIABLE_VALUEdense_20/kernel&variables/8/.ATTRIBUTES/VARIABLE_VALUE
IG
VARIABLE_VALUEdense_20/bias&variables/9/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdense_21/kernel'variables/10/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUEdense_21/bias'variables/11/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdense_22/kernel'variables/12/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUEdense_22/bias'variables/13/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEdense_23/kernel'variables/14/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUEdense_23/bias'variables/15/.ATTRIBUTES/VARIABLE_VALUE

0
1
2

0
1

`0
 
 

!0
"1

!0
"1
 
�
anon_trainable_variables

blayers
cmetrics
dlayer_regularization_losses
elayer_metrics
6	variables
7trainable_variables
8regularization_losses

#0
$1

#0
$1
 
�
fnon_trainable_variables

glayers
hmetrics
ilayer_regularization_losses
jlayer_metrics
:	variables
;trainable_variables
<regularization_losses

%0
&1

%0
&1
 
�
knon_trainable_variables

llayers
mmetrics
nlayer_regularization_losses
olayer_metrics
>	variables
?trainable_variables
@regularization_losses

'0
(1

'0
(1
 
�
pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
B	variables
Ctrainable_variables
Dregularization_losses

)0
*1

)0
*1
 
�
unon_trainable_variables

vlayers
wmetrics
xlayer_regularization_losses
ylayer_metrics
F	variables
Gtrainable_variables
Hregularization_losses

+0
,1

+0
,1
 
�
znon_trainable_variables

{layers
|metrics
}layer_regularization_losses
~layer_metrics
J	variables
Ktrainable_variables
Lregularization_losses

-0
.1

-0
.1
 
�
non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
N	variables
Otrainable_variables
Pregularization_losses
 
1
0
1
2
3
4
5
6
 
 
 

/0
01

/0
01
 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
W	variables
Xtrainable_variables
Yregularization_losses
 

0
 
 
 
8

�total

�count
�	variables
�	keras_api
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
nl
VARIABLE_VALUEAdam/dense_16/kernel/mBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_16/bias/mBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_17/kernel/mBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_17/bias/mBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_18/kernel/mBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_18/bias/mBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_19/kernel/mBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_19/bias/mBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_20/kernel/mBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_20/bias/mBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_21/kernel/mCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_21/bias/mCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_22/kernel/mCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_22/bias/mCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_23/kernel/mCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_23/bias/mCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_16/kernel/vBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_16/bias/vBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_17/kernel/vBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_17/bias/vBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_18/kernel/vBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_18/bias/vBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_19/kernel/vBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_19/bias/vBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
nl
VARIABLE_VALUEAdam/dense_20/kernel/vBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
lj
VARIABLE_VALUEAdam/dense_20/bias/vBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_21/kernel/vCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_21/bias/vCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_22/kernel/vCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_22/bias/vCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/dense_23/kernel/vCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
mk
VARIABLE_VALUEAdam/dense_23/bias/vCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
z
serving_default_input_1Placeholder*'
_output_shapes
:���������*
dtype0*
shape:���������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1dense_16/kerneldense_16/biasdense_17/kerneldense_17/biasdense_18/kerneldense_18/biasdense_19/kerneldense_19/biasdense_20/kerneldense_20/biasdense_21/kerneldense_21/biasdense_22/kerneldense_22/biasVariable
Variable_1
Variable_2dense_23/kerneldense_23/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *-
f(R&
$__inference_signature_wrapper_604570
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameVariable/Read/ReadVariableOpVariable_1/Read/ReadVariableOpVariable_2/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp#dense_16/kernel/Read/ReadVariableOp!dense_16/bias/Read/ReadVariableOp#dense_17/kernel/Read/ReadVariableOp!dense_17/bias/Read/ReadVariableOp#dense_18/kernel/Read/ReadVariableOp!dense_18/bias/Read/ReadVariableOp#dense_19/kernel/Read/ReadVariableOp!dense_19/bias/Read/ReadVariableOp#dense_20/kernel/Read/ReadVariableOp!dense_20/bias/Read/ReadVariableOp#dense_21/kernel/Read/ReadVariableOp!dense_21/bias/Read/ReadVariableOp#dense_22/kernel/Read/ReadVariableOp!dense_22/bias/Read/ReadVariableOp#dense_23/kernel/Read/ReadVariableOp!dense_23/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp*Adam/dense_16/kernel/m/Read/ReadVariableOp(Adam/dense_16/bias/m/Read/ReadVariableOp*Adam/dense_17/kernel/m/Read/ReadVariableOp(Adam/dense_17/bias/m/Read/ReadVariableOp*Adam/dense_18/kernel/m/Read/ReadVariableOp(Adam/dense_18/bias/m/Read/ReadVariableOp*Adam/dense_19/kernel/m/Read/ReadVariableOp(Adam/dense_19/bias/m/Read/ReadVariableOp*Adam/dense_20/kernel/m/Read/ReadVariableOp(Adam/dense_20/bias/m/Read/ReadVariableOp*Adam/dense_21/kernel/m/Read/ReadVariableOp(Adam/dense_21/bias/m/Read/ReadVariableOp*Adam/dense_22/kernel/m/Read/ReadVariableOp(Adam/dense_22/bias/m/Read/ReadVariableOp*Adam/dense_23/kernel/m/Read/ReadVariableOp(Adam/dense_23/bias/m/Read/ReadVariableOp*Adam/dense_16/kernel/v/Read/ReadVariableOp(Adam/dense_16/bias/v/Read/ReadVariableOp*Adam/dense_17/kernel/v/Read/ReadVariableOp(Adam/dense_17/bias/v/Read/ReadVariableOp*Adam/dense_18/kernel/v/Read/ReadVariableOp(Adam/dense_18/bias/v/Read/ReadVariableOp*Adam/dense_19/kernel/v/Read/ReadVariableOp(Adam/dense_19/bias/v/Read/ReadVariableOp*Adam/dense_20/kernel/v/Read/ReadVariableOp(Adam/dense_20/bias/v/Read/ReadVariableOp*Adam/dense_21/kernel/v/Read/ReadVariableOp(Adam/dense_21/bias/v/Read/ReadVariableOp*Adam/dense_22/kernel/v/Read/ReadVariableOp(Adam/dense_22/bias/v/Read/ReadVariableOp*Adam/dense_23/kernel/v/Read/ReadVariableOp(Adam/dense_23/bias/v/Read/ReadVariableOpConst*G
Tin@
>2<	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *(
f#R!
__inference__traced_save_607536
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameVariable
Variable_1
Variable_2	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratedense_16/kerneldense_16/biasdense_17/kerneldense_17/biasdense_18/kerneldense_18/biasdense_19/kerneldense_19/biasdense_20/kerneldense_20/biasdense_21/kerneldense_21/biasdense_22/kerneldense_22/biasdense_23/kerneldense_23/biastotalcountAdam/dense_16/kernel/mAdam/dense_16/bias/mAdam/dense_17/kernel/mAdam/dense_17/bias/mAdam/dense_18/kernel/mAdam/dense_18/bias/mAdam/dense_19/kernel/mAdam/dense_19/bias/mAdam/dense_20/kernel/mAdam/dense_20/bias/mAdam/dense_21/kernel/mAdam/dense_21/bias/mAdam/dense_22/kernel/mAdam/dense_22/bias/mAdam/dense_23/kernel/mAdam/dense_23/bias/mAdam/dense_16/kernel/vAdam/dense_16/bias/vAdam/dense_17/kernel/vAdam/dense_17/bias/vAdam/dense_18/kernel/vAdam/dense_18/bias/vAdam/dense_19/kernel/vAdam/dense_19/bias/vAdam/dense_20/kernel/vAdam/dense_20/bias/vAdam/dense_21/kernel/vAdam/dense_21/bias/vAdam/dense_22/kernel/vAdam/dense_22/bias/vAdam/dense_23/kernel/vAdam/dense_23/bias/v*F
Tin?
=2;*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *+
f&R$
"__inference__traced_restore_607720��4
�0
�
D__inference_dense_17_layer_call_and_return_conditional_losses_601041

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_11_607180D
5dense_21_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOpd
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_21_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_21_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_21/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_15_607339C
5dense_23_bias_regularizer_abs_readvariableop_resource:
identity��,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOpd
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_23_bias_regularizer_abs_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_23_bias_regularizer_abs_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_23/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_4_607040K
7dense_18_kernel_regularizer_abs_readvariableop_resource:
��
identity��.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOpf
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_18_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_18_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_18/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp
�
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493

inputs"
dense_16_600995:	�
dense_16_600997:	�#
dense_17_601042:
��
dense_17_601044:	�#
dense_18_601089:
��
dense_18_601091:	�#
dense_19_601136:
��
dense_19_601138:	�#
dense_20_601183:
��
dense_20_601185:	�#
dense_21_601230:
��
dense_21_601232:	�"
dense_22_601277:	�
dense_22_601279:
identity�� dense_16/StatefulPartitionedCall�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp� dense_17/StatefulPartitionedCall�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp� dense_18/StatefulPartitionedCall�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp� dense_19/StatefulPartitionedCall�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp� dense_20/StatefulPartitionedCall�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp� dense_21/StatefulPartitionedCall�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp� dense_22/StatefulPartitionedCall�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
 dense_16/StatefulPartitionedCallStatefulPartitionedCallinputsdense_16_600995dense_16_600997*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_16_layer_call_and_return_conditional_losses_600994�
 dense_17/StatefulPartitionedCallStatefulPartitionedCall)dense_16/StatefulPartitionedCall:output:0dense_17_601042dense_17_601044*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_17_layer_call_and_return_conditional_losses_601041�
 dense_18/StatefulPartitionedCallStatefulPartitionedCall)dense_17/StatefulPartitionedCall:output:0dense_18_601089dense_18_601091*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_18_layer_call_and_return_conditional_losses_601088�
 dense_19/StatefulPartitionedCallStatefulPartitionedCall)dense_18/StatefulPartitionedCall:output:0dense_19_601136dense_19_601138*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_19_layer_call_and_return_conditional_losses_601135�
 dense_20/StatefulPartitionedCallStatefulPartitionedCall)dense_19/StatefulPartitionedCall:output:0dense_20_601183dense_20_601185*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_20_layer_call_and_return_conditional_losses_601182�
 dense_21/StatefulPartitionedCallStatefulPartitionedCall)dense_20/StatefulPartitionedCall:output:0dense_21_601230dense_21_601232*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_21_layer_call_and_return_conditional_losses_601229�
 dense_22/StatefulPartitionedCallStatefulPartitionedCall)dense_21/StatefulPartitionedCall:output:0dense_22_601277dense_22_601279*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_22_layer_call_and_return_conditional_losses_601276f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_600995*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_600995*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_600997*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_600997*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601042* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601042* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601044*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601044*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601089* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601089* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601091*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601091*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601136* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601136* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601138*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601138*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601183* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601183* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601185*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601185*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601230* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601230* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601232*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601232*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601277*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601277*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601279*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601279*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_22/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_16/StatefulPartitionedCall-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp!^dense_17/StatefulPartitionedCall-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp!^dense_18/StatefulPartitionedCall-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp!^dense_19/StatefulPartitionedCall-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp!^dense_20/StatefulPartitionedCall-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp!^dense_21/StatefulPartitionedCall-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp!^dense_22/StatefulPartitionedCall-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2D
 dense_16/StatefulPartitionedCall dense_16/StatefulPartitionedCall2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2D
 dense_17/StatefulPartitionedCall dense_17/StatefulPartitionedCall2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2D
 dense_18/StatefulPartitionedCall dense_18/StatefulPartitionedCall2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2D
 dense_19/StatefulPartitionedCall dense_19/StatefulPartitionedCall2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2D
 dense_20/StatefulPartitionedCall dense_20/StatefulPartitionedCall2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2D
 dense_21/StatefulPartitionedCall dense_21/StatefulPartitionedCall2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2D
 dense_22/StatefulPartitionedCall dense_22/StatefulPartitionedCall2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_14_607319I
7dense_23_kernel_regularizer_abs_readvariableop_resource:
identity��.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOpf
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_23_kernel_regularizer_abs_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_23_kernel_regularizer_abs_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_23/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_2_607000K
7dense_17_kernel_regularizer_abs_readvariableop_resource:
��
identity��.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOpf
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_17_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_17_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_17/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp
�
�
,__inference_conjugacy_2_layer_call_fn_603575
input_1
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:

unknown_13: 

unknown_14: 

unknown_15: 

unknown_16:

unknown_17:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������: : : *5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603485o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�
�
,__inference_conjugacy_2_layer_call_fn_604616
x
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:

unknown_13: 

unknown_14: 

unknown_15: 

unknown_16:

unknown_17:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������: : : *5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603041o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:J F
'
_output_shapes
:���������

_user_specified_namex
�
�
,__inference_conjugacy_2_layer_call_fn_603085
input_1
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:

unknown_13: 

unknown_14: 

unknown_15: 

unknown_16:

unknown_17:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������: : : *5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603041o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�
�
__inference_loss_fn_12_607200J
7dense_22_kernel_regularizer_abs_readvariableop_resource:	�
identity��.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOpf
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_22_kernel_regularizer_abs_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_22_kernel_regularizer_abs_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_22/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp
�/
�
D__inference_dense_23_layer_call_and_return_conditional_losses_602487

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
-__inference_sequential_4_layer_call_fn_605693

inputs
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_1_606980D
5dense_16_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOpd
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_16_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_16_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_16/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp
�
�
__inference_loss_fn_13_607220C
5dense_22_bias_regularizer_abs_readvariableop_resource:
identity��,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOpd
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_22_bias_regularizer_abs_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_22_bias_regularizer_abs_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_22/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp
��
�
!__inference__wrapped_model_600946
input_1S
@conjugacy_2_sequential_4_dense_16_matmul_readvariableop_resource:	�P
Aconjugacy_2_sequential_4_dense_16_biasadd_readvariableop_resource:	�T
@conjugacy_2_sequential_4_dense_17_matmul_readvariableop_resource:
��P
Aconjugacy_2_sequential_4_dense_17_biasadd_readvariableop_resource:	�T
@conjugacy_2_sequential_4_dense_18_matmul_readvariableop_resource:
��P
Aconjugacy_2_sequential_4_dense_18_biasadd_readvariableop_resource:	�T
@conjugacy_2_sequential_4_dense_19_matmul_readvariableop_resource:
��P
Aconjugacy_2_sequential_4_dense_19_biasadd_readvariableop_resource:	�T
@conjugacy_2_sequential_4_dense_20_matmul_readvariableop_resource:
��P
Aconjugacy_2_sequential_4_dense_20_biasadd_readvariableop_resource:	�T
@conjugacy_2_sequential_4_dense_21_matmul_readvariableop_resource:
��P
Aconjugacy_2_sequential_4_dense_21_biasadd_readvariableop_resource:	�S
@conjugacy_2_sequential_4_dense_22_matmul_readvariableop_resource:	�O
Aconjugacy_2_sequential_4_dense_22_biasadd_readvariableop_resource:-
#conjugacy_2_readvariableop_resource: /
%conjugacy_2_readvariableop_1_resource: /
%conjugacy_2_readvariableop_2_resource: R
@conjugacy_2_sequential_5_dense_23_matmul_readvariableop_resource:O
Aconjugacy_2_sequential_5_dense_23_biasadd_readvariableop_resource:
identity��conjugacy_2/ReadVariableOp�conjugacy_2/ReadVariableOp_1�conjugacy_2/ReadVariableOp_2�8conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOp�9conjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOp�8conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOp�:conjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOp�7conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOp�9conjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOp�
7conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
(conjugacy_2/sequential_4/dense_16/MatMulMatMulinput_1?conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_16/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_16/MatMul:product:0@conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_16/SeluSelu2conjugacy_2/sequential_4/dense_16/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
(conjugacy_2/sequential_4/dense_17/MatMulMatMul4conjugacy_2/sequential_4/dense_16/Selu:activations:0?conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_17/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_17/MatMul:product:0@conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_17/SeluSelu2conjugacy_2/sequential_4/dense_17/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
(conjugacy_2/sequential_4/dense_18/MatMulMatMul4conjugacy_2/sequential_4/dense_17/Selu:activations:0?conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_18/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_18/MatMul:product:0@conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_18/SeluSelu2conjugacy_2/sequential_4/dense_18/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
(conjugacy_2/sequential_4/dense_19/MatMulMatMul4conjugacy_2/sequential_4/dense_18/Selu:activations:0?conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_19/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_19/MatMul:product:0@conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_19/SeluSelu2conjugacy_2/sequential_4/dense_19/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
(conjugacy_2/sequential_4/dense_20/MatMulMatMul4conjugacy_2/sequential_4/dense_19/Selu:activations:0?conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_20/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_20/MatMul:product:0@conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_20/SeluSelu2conjugacy_2/sequential_4/dense_20/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
(conjugacy_2/sequential_4/dense_21/MatMulMatMul4conjugacy_2/sequential_4/dense_20/Selu:activations:0?conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
8conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
)conjugacy_2/sequential_4/dense_21/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_21/MatMul:product:0@conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
&conjugacy_2/sequential_4/dense_21/SeluSelu2conjugacy_2/sequential_4/dense_21/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
7conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
(conjugacy_2/sequential_4/dense_22/MatMulMatMul4conjugacy_2/sequential_4/dense_21/Selu:activations:0?conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
8conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
)conjugacy_2/sequential_4/dense_22/BiasAddBiasAdd2conjugacy_2/sequential_4/dense_22/MatMul:product:0@conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
&conjugacy_2/sequential_4/dense_22/SeluSelu2conjugacy_2/sequential_4/dense_22/BiasAdd:output:0*
T0*'
_output_shapes
:���������p
conjugacy_2/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        r
!conjugacy_2/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       r
!conjugacy_2/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_sliceStridedSlice4conjugacy_2/sequential_4/dense_22/Selu:activations:0(conjugacy_2/strided_slice/stack:output:0*conjugacy_2/strided_slice/stack_1:output:0*conjugacy_2/strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskv
conjugacy_2/ReadVariableOpReadVariableOp#conjugacy_2_readvariableop_resource*
_output_shapes
: *
dtype0�
conjugacy_2/mulMul"conjugacy_2/ReadVariableOp:value:0"conjugacy_2/strided_slice:output:0*
T0*#
_output_shapes
:���������r
!conjugacy_2/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_slice_1StridedSlice4conjugacy_2/sequential_4/dense_22/Selu:activations:0*conjugacy_2/strided_slice_1/stack:output:0,conjugacy_2/strided_slice_1/stack_1:output:0,conjugacy_2/strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskz
conjugacy_2/ReadVariableOp_1ReadVariableOp%conjugacy_2_readvariableop_1_resource*
_output_shapes
: *
dtype0�
conjugacy_2/mul_1Mul$conjugacy_2/ReadVariableOp_1:value:0$conjugacy_2/strided_slice_1:output:0*
T0*#
_output_shapes
:���������r
!conjugacy_2/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_slice_2StridedSlice4conjugacy_2/sequential_4/dense_22/Selu:activations:0*conjugacy_2/strided_slice_2/stack:output:0,conjugacy_2/strided_slice_2/stack_1:output:0,conjugacy_2/strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskz
conjugacy_2/ReadVariableOp_2ReadVariableOp%conjugacy_2_readvariableop_2_resource*
_output_shapes
: *
dtype0�
conjugacy_2/mul_2Mul$conjugacy_2/ReadVariableOp_2:value:0$conjugacy_2/strided_slice_2:output:0*
T0*#
_output_shapes
:����������
conjugacy_2/stackPackconjugacy_2/mul:z:0conjugacy_2/mul_1:z:0conjugacy_2/mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
7conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOpReadVariableOp@conjugacy_2_sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
(conjugacy_2/sequential_5/dense_23/MatMulMatMulconjugacy_2/stack:output:0?conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
8conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOpReadVariableOpAconjugacy_2_sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
)conjugacy_2/sequential_5/dense_23/BiasAddBiasAdd2conjugacy_2/sequential_5/dense_23/MatMul:product:0@conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
9conjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
*conjugacy_2/sequential_5/dense_23/MatMul_1MatMul4conjugacy_2/sequential_4/dense_22/Selu:activations:0Aconjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:conjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
+conjugacy_2/sequential_5/dense_23/BiasAdd_1BiasAdd4conjugacy_2/sequential_5/dense_23/MatMul_1:product:0Bconjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
conjugacy_2/subSubinput_14conjugacy_2/sequential_5/dense_23/BiasAdd_1:output:0*
T0*'
_output_shapes
:���������c
conjugacy_2/SquareSquareconjugacy_2/sub:z:0*
T0*'
_output_shapes
:���������b
conjugacy_2/ConstConst*
_output_shapes
:*
dtype0*
valueB"       m
conjugacy_2/MeanMeanconjugacy_2/Square:y:0conjugacy_2/Const:output:0*
T0*
_output_shapes
: r
!conjugacy_2/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        t
#conjugacy_2/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_slice_3StridedSliceinput_1*conjugacy_2/strided_slice_3/stack:output:0,conjugacy_2/strided_slice_3/stack_1:output:0,conjugacy_2/strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskX
conjugacy_2/mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
conjugacy_2/mul_3Mulconjugacy_2/mul_3/x:output:0$conjugacy_2/strided_slice_3:output:0*
T0*#
_output_shapes
:���������r
!conjugacy_2/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_slice_4StridedSliceinput_1*conjugacy_2/strided_slice_4/stack:output:0,conjugacy_2/strided_slice_4/stack_1:output:0,conjugacy_2/strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskr
!conjugacy_2/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        t
#conjugacy_2/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       t
#conjugacy_2/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
conjugacy_2/strided_slice_5StridedSliceinput_1*conjugacy_2/strided_slice_5/stack:output:0,conjugacy_2/strided_slice_5/stack_1:output:0,conjugacy_2/strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskr
conjugacy_2/Square_1Square$conjugacy_2/strided_slice_5:output:0*
T0*#
_output_shapes
:����������
conjugacy_2/sub_1Sub$conjugacy_2/strided_slice_4:output:0conjugacy_2/Square_1:y:0*
T0*#
_output_shapes
:���������X
conjugacy_2/mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?{
conjugacy_2/mul_4Mulconjugacy_2/mul_4/x:output:0conjugacy_2/sub_1:z:0*
T0*#
_output_shapes
:����������
conjugacy_2/stack_1Packconjugacy_2/mul_3:z:0conjugacy_2/mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
conjugacy_2/sub_2Sub2conjugacy_2/sequential_5/dense_23/BiasAdd:output:0conjugacy_2/stack_1:output:0*
T0*'
_output_shapes
:���������g
conjugacy_2/Square_2Squareconjugacy_2/sub_2:z:0*
T0*'
_output_shapes
:���������d
conjugacy_2/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       s
conjugacy_2/Mean_1Meanconjugacy_2/Square_2:y:0conjugacy_2/Const_1:output:0*
T0*
_output_shapes
: �
9conjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
*conjugacy_2/sequential_4/dense_16/MatMul_1MatMulconjugacy_2/stack_1:output:0Aconjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_16/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_16/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_16/Selu_1Selu4conjugacy_2/sequential_4/dense_16/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
*conjugacy_2/sequential_4/dense_17/MatMul_1MatMul6conjugacy_2/sequential_4/dense_16/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_17/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_17/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_17/Selu_1Selu4conjugacy_2/sequential_4/dense_17/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
*conjugacy_2/sequential_4/dense_18/MatMul_1MatMul6conjugacy_2/sequential_4/dense_17/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_18/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_18/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_18/Selu_1Selu4conjugacy_2/sequential_4/dense_18/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
*conjugacy_2/sequential_4/dense_19/MatMul_1MatMul6conjugacy_2/sequential_4/dense_18/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_19/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_19/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_19/Selu_1Selu4conjugacy_2/sequential_4/dense_19/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
*conjugacy_2/sequential_4/dense_20/MatMul_1MatMul6conjugacy_2/sequential_4/dense_19/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_20/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_20/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_20/Selu_1Selu4conjugacy_2/sequential_4/dense_20/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
*conjugacy_2/sequential_4/dense_21/MatMul_1MatMul6conjugacy_2/sequential_4/dense_20/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
:conjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
+conjugacy_2/sequential_4/dense_21/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_21/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
(conjugacy_2/sequential_4/dense_21/Selu_1Selu4conjugacy_2/sequential_4/dense_21/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
9conjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOpReadVariableOp@conjugacy_2_sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
*conjugacy_2/sequential_4/dense_22/MatMul_1MatMul6conjugacy_2/sequential_4/dense_21/Selu_1:activations:0Aconjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
:conjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOpReadVariableOpAconjugacy_2_sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
+conjugacy_2/sequential_4/dense_22/BiasAdd_1BiasAdd4conjugacy_2/sequential_4/dense_22/MatMul_1:product:0Bconjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
(conjugacy_2/sequential_4/dense_22/Selu_1Selu4conjugacy_2/sequential_4/dense_22/BiasAdd_1:output:0*
T0*'
_output_shapes
:����������
conjugacy_2/sub_3Subconjugacy_2/stack:output:06conjugacy_2/sequential_4/dense_22/Selu_1:activations:0*
T0*'
_output_shapes
:���������g
conjugacy_2/Square_3Squareconjugacy_2/sub_3:z:0*
T0*'
_output_shapes
:���������d
conjugacy_2/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       s
conjugacy_2/Mean_2Meanconjugacy_2/Square_3:y:0conjugacy_2/Const_2:output:0*
T0*
_output_shapes
: �
IdentityIdentity2conjugacy_2/sequential_5/dense_23/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^conjugacy_2/ReadVariableOp^conjugacy_2/ReadVariableOp_1^conjugacy_2/ReadVariableOp_29^conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOp:^conjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOp9^conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOp;^conjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOp8^conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOp:^conjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 28
conjugacy_2/ReadVariableOpconjugacy_2/ReadVariableOp2<
conjugacy_2/ReadVariableOp_1conjugacy_2/ReadVariableOp_12<
conjugacy_2/ReadVariableOp_2conjugacy_2/ReadVariableOp_22t
8conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_16/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_16/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_16/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_16/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_17/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_17/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_17/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_17/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_18/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_18/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_18/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_18/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_19/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_19/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_19/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_19/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_20/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_20/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_20/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_20/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_21/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_21/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_21/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_21/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOp8conjugacy_2/sequential_4/dense_22/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_4/dense_22/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOp7conjugacy_2/sequential_4/dense_22/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOp9conjugacy_2/sequential_4/dense_22/MatMul_1/ReadVariableOp2t
8conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOp8conjugacy_2/sequential_5/dense_23/BiasAdd/ReadVariableOp2x
:conjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOp:conjugacy_2/sequential_5/dense_23/BiasAdd_1/ReadVariableOp2r
7conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOp7conjugacy_2/sequential_5/dense_23/MatMul/ReadVariableOp2v
9conjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOp9conjugacy_2/sequential_5/dense_23/MatMul_1/ReadVariableOp:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�0
�
D__inference_dense_21_layer_call_and_return_conditional_losses_601229

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�1
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_606380

inputs9
'dense_23_matmul_readvariableop_resource:6
(dense_23_biasadd_readvariableop_resource:
identity��dense_23/BiasAdd/ReadVariableOp�dense_23/MatMul/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
dense_23/MatMul/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0{
dense_23/MatMulMatMulinputs&dense_23/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
dense_23/BiasAdd/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_23/BiasAddBiasAdddense_23/MatMul:product:0'dense_23/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: h
IdentityIdentitydense_23/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp ^dense_23/BiasAdd/ReadVariableOp^dense_23/MatMul/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2B
dense_23/BiasAdd/ReadVariableOpdense_23/BiasAdd/ReadVariableOp2@
dense_23/MatMul/ReadVariableOpdense_23/MatMul/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�0
�
D__inference_dense_20_layer_call_and_return_conditional_losses_606780

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603927
input_1&
sequential_4_603578:	�"
sequential_4_603580:	�'
sequential_4_603582:
��"
sequential_4_603584:	�'
sequential_4_603586:
��"
sequential_4_603588:	�'
sequential_4_603590:
��"
sequential_4_603592:	�'
sequential_4_603594:
��"
sequential_4_603596:	�'
sequential_4_603598:
��"
sequential_4_603600:	�&
sequential_4_603602:	�!
sequential_4_603604:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: %
sequential_5_603629:!
sequential_5_603631:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�&sequential_4/StatefulPartitionedCall_1�$sequential_5/StatefulPartitionedCall�&sequential_5/StatefulPartitionedCall_1�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallinput_1sequential_4_603578sequential_4_603580sequential_4_603582sequential_4_603584sequential_4_603586sequential_4_603588sequential_4_603590sequential_4_603592sequential_4_603594sequential_4_603596sequential_4_603598sequential_4_603600sequential_4_603602sequential_4_603604*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
$sequential_5/StatefulPartitionedCallStatefulPartitionedCallstack:output:0sequential_5_603629sequential_5_603631*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524�
&sequential_5/StatefulPartitionedCall_1StatefulPartitionedCall-sequential_4/StatefulPartitionedCall:output:0sequential_5_603629sequential_5_603631*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524v
subSubinput_1/sequential_5/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSliceinput_1strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSliceinput_1strided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSliceinput_1strided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������
sub_2Sub-sequential_5/StatefulPartitionedCall:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCallstack_1:output:0sequential_4_603578sequential_4_603580sequential_4_603582sequential_4_603584sequential_4_603586sequential_4_603588sequential_4_603590sequential_4_603592sequential_4_603594sequential_4_603596sequential_4_603598sequential_4_603600sequential_4_603602sequential_4_603604*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493
sub_3Substack:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603578*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603578*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603580*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603580*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603582* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603582* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603584*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603584*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603586* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603586* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603588*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603588*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603590* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603590* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603592*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603592*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603594* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603594* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603596*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603596*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603598* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603598* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603600*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603600*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603602*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603602*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603604*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603604*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603629*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603629*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603631*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603631*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: |
IdentityIdentity-sequential_5/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1%^sequential_5/StatefulPartitionedCall'^sequential_5/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_12L
$sequential_5/StatefulPartitionedCall$sequential_5/StatefulPartitionedCall2P
&sequential_5/StatefulPartitionedCall_1&sequential_5/StatefulPartitionedCall_1:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
܍
�#
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605056
xG
4sequential_4_dense_16_matmul_readvariableop_resource:	�D
5sequential_4_dense_16_biasadd_readvariableop_resource:	�H
4sequential_4_dense_17_matmul_readvariableop_resource:
��D
5sequential_4_dense_17_biasadd_readvariableop_resource:	�H
4sequential_4_dense_18_matmul_readvariableop_resource:
��D
5sequential_4_dense_18_biasadd_readvariableop_resource:	�H
4sequential_4_dense_19_matmul_readvariableop_resource:
��D
5sequential_4_dense_19_biasadd_readvariableop_resource:	�H
4sequential_4_dense_20_matmul_readvariableop_resource:
��D
5sequential_4_dense_20_biasadd_readvariableop_resource:	�H
4sequential_4_dense_21_matmul_readvariableop_resource:
��D
5sequential_4_dense_21_biasadd_readvariableop_resource:	�G
4sequential_4_dense_22_matmul_readvariableop_resource:	�C
5sequential_4_dense_22_biasadd_readvariableop_resource:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: F
4sequential_5_dense_23_matmul_readvariableop_resource:C
5sequential_5_dense_23_biasadd_readvariableop_resource:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�,sequential_4/dense_16/BiasAdd/ReadVariableOp�.sequential_4/dense_16/BiasAdd_1/ReadVariableOp�+sequential_4/dense_16/MatMul/ReadVariableOp�-sequential_4/dense_16/MatMul_1/ReadVariableOp�,sequential_4/dense_17/BiasAdd/ReadVariableOp�.sequential_4/dense_17/BiasAdd_1/ReadVariableOp�+sequential_4/dense_17/MatMul/ReadVariableOp�-sequential_4/dense_17/MatMul_1/ReadVariableOp�,sequential_4/dense_18/BiasAdd/ReadVariableOp�.sequential_4/dense_18/BiasAdd_1/ReadVariableOp�+sequential_4/dense_18/MatMul/ReadVariableOp�-sequential_4/dense_18/MatMul_1/ReadVariableOp�,sequential_4/dense_19/BiasAdd/ReadVariableOp�.sequential_4/dense_19/BiasAdd_1/ReadVariableOp�+sequential_4/dense_19/MatMul/ReadVariableOp�-sequential_4/dense_19/MatMul_1/ReadVariableOp�,sequential_4/dense_20/BiasAdd/ReadVariableOp�.sequential_4/dense_20/BiasAdd_1/ReadVariableOp�+sequential_4/dense_20/MatMul/ReadVariableOp�-sequential_4/dense_20/MatMul_1/ReadVariableOp�,sequential_4/dense_21/BiasAdd/ReadVariableOp�.sequential_4/dense_21/BiasAdd_1/ReadVariableOp�+sequential_4/dense_21/MatMul/ReadVariableOp�-sequential_4/dense_21/MatMul_1/ReadVariableOp�,sequential_4/dense_22/BiasAdd/ReadVariableOp�.sequential_4/dense_22/BiasAdd_1/ReadVariableOp�+sequential_4/dense_22/MatMul/ReadVariableOp�-sequential_4/dense_22/MatMul_1/ReadVariableOp�,sequential_5/dense_23/BiasAdd/ReadVariableOp�.sequential_5/dense_23/BiasAdd_1/ReadVariableOp�+sequential_5/dense_23/MatMul/ReadVariableOp�-sequential_5/dense_23/MatMul_1/ReadVariableOp�
+sequential_4/dense_16/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_16/MatMulMatMulx3sequential_4/dense_16/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_16/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_16/BiasAddBiasAdd&sequential_4/dense_16/MatMul:product:04sequential_4/dense_16/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_16/SeluSelu&sequential_4/dense_16/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_17/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_17/MatMulMatMul(sequential_4/dense_16/Selu:activations:03sequential_4/dense_17/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_17/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_17/BiasAddBiasAdd&sequential_4/dense_17/MatMul:product:04sequential_4/dense_17/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_17/SeluSelu&sequential_4/dense_17/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_18/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_18/MatMulMatMul(sequential_4/dense_17/Selu:activations:03sequential_4/dense_18/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_18/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_18/BiasAddBiasAdd&sequential_4/dense_18/MatMul:product:04sequential_4/dense_18/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_18/SeluSelu&sequential_4/dense_18/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_19/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_19/MatMulMatMul(sequential_4/dense_18/Selu:activations:03sequential_4/dense_19/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_19/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_19/BiasAddBiasAdd&sequential_4/dense_19/MatMul:product:04sequential_4/dense_19/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_19/SeluSelu&sequential_4/dense_19/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_20/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_20/MatMulMatMul(sequential_4/dense_19/Selu:activations:03sequential_4/dense_20/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_20/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_20/BiasAddBiasAdd&sequential_4/dense_20/MatMul:product:04sequential_4/dense_20/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_20/SeluSelu&sequential_4/dense_20/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_21/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_21/MatMulMatMul(sequential_4/dense_20/Selu:activations:03sequential_4/dense_21/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_21/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_21/BiasAddBiasAdd&sequential_4/dense_21/MatMul:product:04sequential_4/dense_21/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_21/SeluSelu&sequential_4/dense_21/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_22/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_22/MatMulMatMul(sequential_4/dense_21/Selu:activations:03sequential_4/dense_22/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
,sequential_4/dense_22/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_4/dense_22/BiasAddBiasAdd&sequential_4/dense_22/MatMul:product:04sequential_4/dense_22/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������|
sequential_4/dense_22/SeluSelu&sequential_4/dense_22/BiasAdd:output:0*
T0*'
_output_shapes
:���������d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
+sequential_5/dense_23/MatMul/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sequential_5/dense_23/MatMulMatMulstack:output:03sequential_5/dense_23/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
,sequential_5/dense_23/BiasAdd/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_5/dense_23/BiasAddBiasAdd&sequential_5/dense_23/MatMul:product:04sequential_5/dense_23/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
-sequential_5/dense_23/MatMul_1/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sequential_5/dense_23/MatMul_1MatMul(sequential_4/dense_22/Selu:activations:05sequential_5/dense_23/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_5/dense_23/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_5/dense_23/BiasAdd_1BiasAdd(sequential_5/dense_23/MatMul_1:product:06sequential_5/dense_23/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������i
subSubx(sequential_5/dense_23/BiasAdd_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSlicexstrided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSlicexstrided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSlicexstrided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������x
sub_2Sub&sequential_5/dense_23/BiasAdd:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
-sequential_4/dense_16/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_16/MatMul_1MatMulstack_1:output:05sequential_4/dense_16/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_16/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_16/BiasAdd_1BiasAdd(sequential_4/dense_16/MatMul_1:product:06sequential_4/dense_16/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_16/Selu_1Selu(sequential_4/dense_16/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_17/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_17/MatMul_1MatMul*sequential_4/dense_16/Selu_1:activations:05sequential_4/dense_17/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_17/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_17/BiasAdd_1BiasAdd(sequential_4/dense_17/MatMul_1:product:06sequential_4/dense_17/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_17/Selu_1Selu(sequential_4/dense_17/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_18/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_18/MatMul_1MatMul*sequential_4/dense_17/Selu_1:activations:05sequential_4/dense_18/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_18/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_18/BiasAdd_1BiasAdd(sequential_4/dense_18/MatMul_1:product:06sequential_4/dense_18/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_18/Selu_1Selu(sequential_4/dense_18/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_19/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_19/MatMul_1MatMul*sequential_4/dense_18/Selu_1:activations:05sequential_4/dense_19/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_19/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_19/BiasAdd_1BiasAdd(sequential_4/dense_19/MatMul_1:product:06sequential_4/dense_19/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_19/Selu_1Selu(sequential_4/dense_19/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_20/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_20/MatMul_1MatMul*sequential_4/dense_19/Selu_1:activations:05sequential_4/dense_20/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_20/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_20/BiasAdd_1BiasAdd(sequential_4/dense_20/MatMul_1:product:06sequential_4/dense_20/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_20/Selu_1Selu(sequential_4/dense_20/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_21/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_21/MatMul_1MatMul*sequential_4/dense_20/Selu_1:activations:05sequential_4/dense_21/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_21/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_21/BiasAdd_1BiasAdd(sequential_4/dense_21/MatMul_1:product:06sequential_4/dense_21/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_21/Selu_1Selu(sequential_4/dense_21/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_22/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_22/MatMul_1MatMul*sequential_4/dense_21/Selu_1:activations:05sequential_4/dense_22/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_4/dense_22/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_4/dense_22/BiasAdd_1BiasAdd(sequential_4/dense_22/MatMul_1:product:06sequential_4/dense_22/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
sequential_4/dense_22/Selu_1Selu(sequential_4/dense_22/BiasAdd_1:output:0*
T0*'
_output_shapes
:���������z
sub_3Substack:output:0*sequential_4/dense_22/Selu_1:activations:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: u
IdentityIdentity&sequential_5/dense_23/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp-^sequential_4/dense_16/BiasAdd/ReadVariableOp/^sequential_4/dense_16/BiasAdd_1/ReadVariableOp,^sequential_4/dense_16/MatMul/ReadVariableOp.^sequential_4/dense_16/MatMul_1/ReadVariableOp-^sequential_4/dense_17/BiasAdd/ReadVariableOp/^sequential_4/dense_17/BiasAdd_1/ReadVariableOp,^sequential_4/dense_17/MatMul/ReadVariableOp.^sequential_4/dense_17/MatMul_1/ReadVariableOp-^sequential_4/dense_18/BiasAdd/ReadVariableOp/^sequential_4/dense_18/BiasAdd_1/ReadVariableOp,^sequential_4/dense_18/MatMul/ReadVariableOp.^sequential_4/dense_18/MatMul_1/ReadVariableOp-^sequential_4/dense_19/BiasAdd/ReadVariableOp/^sequential_4/dense_19/BiasAdd_1/ReadVariableOp,^sequential_4/dense_19/MatMul/ReadVariableOp.^sequential_4/dense_19/MatMul_1/ReadVariableOp-^sequential_4/dense_20/BiasAdd/ReadVariableOp/^sequential_4/dense_20/BiasAdd_1/ReadVariableOp,^sequential_4/dense_20/MatMul/ReadVariableOp.^sequential_4/dense_20/MatMul_1/ReadVariableOp-^sequential_4/dense_21/BiasAdd/ReadVariableOp/^sequential_4/dense_21/BiasAdd_1/ReadVariableOp,^sequential_4/dense_21/MatMul/ReadVariableOp.^sequential_4/dense_21/MatMul_1/ReadVariableOp-^sequential_4/dense_22/BiasAdd/ReadVariableOp/^sequential_4/dense_22/BiasAdd_1/ReadVariableOp,^sequential_4/dense_22/MatMul/ReadVariableOp.^sequential_4/dense_22/MatMul_1/ReadVariableOp-^sequential_5/dense_23/BiasAdd/ReadVariableOp/^sequential_5/dense_23/BiasAdd_1/ReadVariableOp,^sequential_5/dense_23/MatMul/ReadVariableOp.^sequential_5/dense_23/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2\
,sequential_4/dense_16/BiasAdd/ReadVariableOp,sequential_4/dense_16/BiasAdd/ReadVariableOp2`
.sequential_4/dense_16/BiasAdd_1/ReadVariableOp.sequential_4/dense_16/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_16/MatMul/ReadVariableOp+sequential_4/dense_16/MatMul/ReadVariableOp2^
-sequential_4/dense_16/MatMul_1/ReadVariableOp-sequential_4/dense_16/MatMul_1/ReadVariableOp2\
,sequential_4/dense_17/BiasAdd/ReadVariableOp,sequential_4/dense_17/BiasAdd/ReadVariableOp2`
.sequential_4/dense_17/BiasAdd_1/ReadVariableOp.sequential_4/dense_17/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_17/MatMul/ReadVariableOp+sequential_4/dense_17/MatMul/ReadVariableOp2^
-sequential_4/dense_17/MatMul_1/ReadVariableOp-sequential_4/dense_17/MatMul_1/ReadVariableOp2\
,sequential_4/dense_18/BiasAdd/ReadVariableOp,sequential_4/dense_18/BiasAdd/ReadVariableOp2`
.sequential_4/dense_18/BiasAdd_1/ReadVariableOp.sequential_4/dense_18/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_18/MatMul/ReadVariableOp+sequential_4/dense_18/MatMul/ReadVariableOp2^
-sequential_4/dense_18/MatMul_1/ReadVariableOp-sequential_4/dense_18/MatMul_1/ReadVariableOp2\
,sequential_4/dense_19/BiasAdd/ReadVariableOp,sequential_4/dense_19/BiasAdd/ReadVariableOp2`
.sequential_4/dense_19/BiasAdd_1/ReadVariableOp.sequential_4/dense_19/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_19/MatMul/ReadVariableOp+sequential_4/dense_19/MatMul/ReadVariableOp2^
-sequential_4/dense_19/MatMul_1/ReadVariableOp-sequential_4/dense_19/MatMul_1/ReadVariableOp2\
,sequential_4/dense_20/BiasAdd/ReadVariableOp,sequential_4/dense_20/BiasAdd/ReadVariableOp2`
.sequential_4/dense_20/BiasAdd_1/ReadVariableOp.sequential_4/dense_20/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_20/MatMul/ReadVariableOp+sequential_4/dense_20/MatMul/ReadVariableOp2^
-sequential_4/dense_20/MatMul_1/ReadVariableOp-sequential_4/dense_20/MatMul_1/ReadVariableOp2\
,sequential_4/dense_21/BiasAdd/ReadVariableOp,sequential_4/dense_21/BiasAdd/ReadVariableOp2`
.sequential_4/dense_21/BiasAdd_1/ReadVariableOp.sequential_4/dense_21/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_21/MatMul/ReadVariableOp+sequential_4/dense_21/MatMul/ReadVariableOp2^
-sequential_4/dense_21/MatMul_1/ReadVariableOp-sequential_4/dense_21/MatMul_1/ReadVariableOp2\
,sequential_4/dense_22/BiasAdd/ReadVariableOp,sequential_4/dense_22/BiasAdd/ReadVariableOp2`
.sequential_4/dense_22/BiasAdd_1/ReadVariableOp.sequential_4/dense_22/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_22/MatMul/ReadVariableOp+sequential_4/dense_22/MatMul/ReadVariableOp2^
-sequential_4/dense_22/MatMul_1/ReadVariableOp-sequential_4/dense_22/MatMul_1/ReadVariableOp2\
,sequential_5/dense_23/BiasAdd/ReadVariableOp,sequential_5/dense_23/BiasAdd/ReadVariableOp2`
.sequential_5/dense_23/BiasAdd_1/ReadVariableOp.sequential_5/dense_23/BiasAdd_1/ReadVariableOp2Z
+sequential_5/dense_23/MatMul/ReadVariableOp+sequential_5/dense_23/MatMul/ReadVariableOp2^
-sequential_5/dense_23/MatMul_1/ReadVariableOp-sequential_5/dense_23/MatMul_1/ReadVariableOp:J F
'
_output_shapes
:���������

_user_specified_namex
݁
�
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603485
x&
sequential_4_603136:	�"
sequential_4_603138:	�'
sequential_4_603140:
��"
sequential_4_603142:	�'
sequential_4_603144:
��"
sequential_4_603146:	�'
sequential_4_603148:
��"
sequential_4_603150:	�'
sequential_4_603152:
��"
sequential_4_603154:	�'
sequential_4_603156:
��"
sequential_4_603158:	�&
sequential_4_603160:	�!
sequential_4_603162:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: %
sequential_5_603187:!
sequential_5_603189:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�&sequential_4/StatefulPartitionedCall_1�$sequential_5/StatefulPartitionedCall�&sequential_5/StatefulPartitionedCall_1�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallxsequential_4_603136sequential_4_603138sequential_4_603140sequential_4_603142sequential_4_603144sequential_4_603146sequential_4_603148sequential_4_603150sequential_4_603152sequential_4_603154sequential_4_603156sequential_4_603158sequential_4_603160sequential_4_603162*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
$sequential_5/StatefulPartitionedCallStatefulPartitionedCallstack:output:0sequential_5_603187sequential_5_603189*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591�
&sequential_5/StatefulPartitionedCall_1StatefulPartitionedCall-sequential_4/StatefulPartitionedCall:output:0sequential_5_603187sequential_5_603189*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591p
subSubx/sequential_5/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSlicexstrided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSlicexstrided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSlicexstrided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������
sub_2Sub-sequential_5/StatefulPartitionedCall:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCallstack_1:output:0sequential_4_603136sequential_4_603138sequential_4_603140sequential_4_603142sequential_4_603144sequential_4_603146sequential_4_603148sequential_4_603150sequential_4_603152sequential_4_603154sequential_4_603156sequential_4_603158sequential_4_603160sequential_4_603162*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878
sub_3Substack:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603136*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603136*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603138*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603138*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603140* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603140* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603142*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603142*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603144* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603144* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603146*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603146*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603148* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603148* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603150*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603150*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603152* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603152* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603154*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603154*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603156* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603156* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603158*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603158*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603160*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603160*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603162*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603162*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603187*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603187*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603189*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603189*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: |
IdentityIdentity-sequential_5/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1%^sequential_5/StatefulPartitionedCall'^sequential_5/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_12L
$sequential_5/StatefulPartitionedCall$sequential_5/StatefulPartitionedCall2P
&sequential_5/StatefulPartitionedCall_1&sequential_5/StatefulPartitionedCall_1:J F
'
_output_shapes
:���������

_user_specified_namex
�
�
)__inference_dense_16_layer_call_fn_606419

inputs
unknown:	�
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_16_layer_call_and_return_conditional_losses_600994p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�0
�
D__inference_dense_19_layer_call_and_return_conditional_losses_606700

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
-__inference_sequential_5_layer_call_fn_602531
dense_23_input
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_23_inputunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:���������
(
_user_specified_namedense_23_input
�0
�
D__inference_dense_16_layer_call_and_return_conditional_losses_606460

inputs1
matmul_readvariableop_resource:	�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�1
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_606340

inputs9
'dense_23_matmul_readvariableop_resource:6
(dense_23_biasadd_readvariableop_resource:
identity��dense_23/BiasAdd/ReadVariableOp�dense_23/MatMul/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
dense_23/MatMul/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0{
dense_23/MatMulMatMulinputs&dense_23/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
dense_23/BiasAdd/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_23/BiasAddBiasAdddense_23/MatMul:product:0'dense_23/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: h
IdentityIdentitydense_23/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp ^dense_23/BiasAdd/ReadVariableOp^dense_23/MatMul/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2B
dense_23/BiasAdd/ReadVariableOpdense_23/BiasAdd/ReadVariableOp2@
dense_23/MatMul/ReadVariableOpdense_23/MatMul/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
)__inference_dense_21_layer_call_fn_606819

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_21_layer_call_and_return_conditional_losses_601229p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
��
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_602440
dense_16_input"
dense_16_602194:	�
dense_16_602196:	�#
dense_17_602199:
��
dense_17_602201:	�#
dense_18_602204:
��
dense_18_602206:	�#
dense_19_602209:
��
dense_19_602211:	�#
dense_20_602214:
��
dense_20_602216:	�#
dense_21_602219:
��
dense_21_602221:	�"
dense_22_602224:	�
dense_22_602226:
identity�� dense_16/StatefulPartitionedCall�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp� dense_17/StatefulPartitionedCall�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp� dense_18/StatefulPartitionedCall�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp� dense_19/StatefulPartitionedCall�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp� dense_20/StatefulPartitionedCall�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp� dense_21/StatefulPartitionedCall�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp� dense_22/StatefulPartitionedCall�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
 dense_16/StatefulPartitionedCallStatefulPartitionedCalldense_16_inputdense_16_602194dense_16_602196*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_16_layer_call_and_return_conditional_losses_600994�
 dense_17/StatefulPartitionedCallStatefulPartitionedCall)dense_16/StatefulPartitionedCall:output:0dense_17_602199dense_17_602201*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_17_layer_call_and_return_conditional_losses_601041�
 dense_18/StatefulPartitionedCallStatefulPartitionedCall)dense_17/StatefulPartitionedCall:output:0dense_18_602204dense_18_602206*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_18_layer_call_and_return_conditional_losses_601088�
 dense_19/StatefulPartitionedCallStatefulPartitionedCall)dense_18/StatefulPartitionedCall:output:0dense_19_602209dense_19_602211*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_19_layer_call_and_return_conditional_losses_601135�
 dense_20/StatefulPartitionedCallStatefulPartitionedCall)dense_19/StatefulPartitionedCall:output:0dense_20_602214dense_20_602216*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_20_layer_call_and_return_conditional_losses_601182�
 dense_21/StatefulPartitionedCallStatefulPartitionedCall)dense_20/StatefulPartitionedCall:output:0dense_21_602219dense_21_602221*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_21_layer_call_and_return_conditional_losses_601229�
 dense_22/StatefulPartitionedCallStatefulPartitionedCall)dense_21/StatefulPartitionedCall:output:0dense_22_602224dense_22_602226*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_22_layer_call_and_return_conditional_losses_601276f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_602194*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_602194*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_602196*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_602196*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_602199* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_602199* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_602201*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_602201*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_602204* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_602204* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_602206*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_602206*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_602209* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_602209* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_602211*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_602211*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_602214* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_602214* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_602216*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_602216*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_602219* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_602219* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_602221*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_602221*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_602224*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_602224*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_602226*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_602226*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_22/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_16/StatefulPartitionedCall-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp!^dense_17/StatefulPartitionedCall-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp!^dense_18/StatefulPartitionedCall-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp!^dense_19/StatefulPartitionedCall-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp!^dense_20/StatefulPartitionedCall-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp!^dense_21/StatefulPartitionedCall-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp!^dense_22/StatefulPartitionedCall-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2D
 dense_16/StatefulPartitionedCall dense_16/StatefulPartitionedCall2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2D
 dense_17/StatefulPartitionedCall dense_17/StatefulPartitionedCall2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2D
 dense_18/StatefulPartitionedCall dense_18/StatefulPartitionedCall2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2D
 dense_19/StatefulPartitionedCall dense_19/StatefulPartitionedCall2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2D
 dense_20/StatefulPartitionedCall dense_20/StatefulPartitionedCall2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2D
 dense_21/StatefulPartitionedCall dense_21/StatefulPartitionedCall2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2D
 dense_22/StatefulPartitionedCall dense_22/StatefulPartitionedCall2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:W S
'
_output_shapes
:���������
(
_user_specified_namedense_16_input
�
�
__inference_loss_fn_7_607100D
5dense_19_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOpd
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_19_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_19_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_19/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp
�
�
-__inference_sequential_4_layer_call_fn_601524
dense_16_input
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_16_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:���������
(
_user_specified_namedense_16_input
��
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_602191
dense_16_input"
dense_16_601945:	�
dense_16_601947:	�#
dense_17_601950:
��
dense_17_601952:	�#
dense_18_601955:
��
dense_18_601957:	�#
dense_19_601960:
��
dense_19_601962:	�#
dense_20_601965:
��
dense_20_601967:	�#
dense_21_601970:
��
dense_21_601972:	�"
dense_22_601975:	�
dense_22_601977:
identity�� dense_16/StatefulPartitionedCall�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp� dense_17/StatefulPartitionedCall�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp� dense_18/StatefulPartitionedCall�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp� dense_19/StatefulPartitionedCall�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp� dense_20/StatefulPartitionedCall�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp� dense_21/StatefulPartitionedCall�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp� dense_22/StatefulPartitionedCall�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
 dense_16/StatefulPartitionedCallStatefulPartitionedCalldense_16_inputdense_16_601945dense_16_601947*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_16_layer_call_and_return_conditional_losses_600994�
 dense_17/StatefulPartitionedCallStatefulPartitionedCall)dense_16/StatefulPartitionedCall:output:0dense_17_601950dense_17_601952*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_17_layer_call_and_return_conditional_losses_601041�
 dense_18/StatefulPartitionedCallStatefulPartitionedCall)dense_17/StatefulPartitionedCall:output:0dense_18_601955dense_18_601957*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_18_layer_call_and_return_conditional_losses_601088�
 dense_19/StatefulPartitionedCallStatefulPartitionedCall)dense_18/StatefulPartitionedCall:output:0dense_19_601960dense_19_601962*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_19_layer_call_and_return_conditional_losses_601135�
 dense_20/StatefulPartitionedCallStatefulPartitionedCall)dense_19/StatefulPartitionedCall:output:0dense_20_601965dense_20_601967*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_20_layer_call_and_return_conditional_losses_601182�
 dense_21/StatefulPartitionedCallStatefulPartitionedCall)dense_20/StatefulPartitionedCall:output:0dense_21_601970dense_21_601972*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_21_layer_call_and_return_conditional_losses_601229�
 dense_22/StatefulPartitionedCallStatefulPartitionedCall)dense_21/StatefulPartitionedCall:output:0dense_22_601975dense_22_601977*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_22_layer_call_and_return_conditional_losses_601276f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_601945*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_601945*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_601947*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_601947*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601950* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601950* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601952*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601952*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601955* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601955* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601957*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601957*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601960* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601960* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601962*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601962*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601965* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601965* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601967*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601967*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601970* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601970* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601972*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601972*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601975*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601975*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601977*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601977*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_22/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_16/StatefulPartitionedCall-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp!^dense_17/StatefulPartitionedCall-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp!^dense_18/StatefulPartitionedCall-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp!^dense_19/StatefulPartitionedCall-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp!^dense_20/StatefulPartitionedCall-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp!^dense_21/StatefulPartitionedCall-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp!^dense_22/StatefulPartitionedCall-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2D
 dense_16/StatefulPartitionedCall dense_16/StatefulPartitionedCall2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2D
 dense_17/StatefulPartitionedCall dense_17/StatefulPartitionedCall2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2D
 dense_18/StatefulPartitionedCall dense_18/StatefulPartitionedCall2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2D
 dense_19/StatefulPartitionedCall dense_19/StatefulPartitionedCall2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2D
 dense_20/StatefulPartitionedCall dense_20/StatefulPartitionedCall2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2D
 dense_21/StatefulPartitionedCall dense_21/StatefulPartitionedCall2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2D
 dense_22/StatefulPartitionedCall dense_22/StatefulPartitionedCall2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:W S
'
_output_shapes
:���������
(
_user_specified_namedense_16_input
�
�
-__inference_sequential_5_layer_call_fn_602607
dense_23_input
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_23_inputunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:���������
(
_user_specified_namedense_23_input
�
�
-__inference_sequential_5_layer_call_fn_606291

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
,__inference_conjugacy_2_layer_call_fn_604662
x
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:

unknown_13: 

unknown_14: 

unknown_15: 

unknown_16:

unknown_17:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallxunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17*
Tin
2*
Tout
2*
_collective_manager_ids
 *-
_output_shapes
:���������: : : *5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603485o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:J F
'
_output_shapes
:���������

_user_specified_namex
�
�
__inference_loss_fn_8_607120K
7dense_20_kernel_regularizer_abs_readvariableop_resource:
��
identity��.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOpf
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_20_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_20_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_20/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp
�.
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_602646
dense_23_input!
dense_23_602610:
dense_23_602612:
identity�� dense_23/StatefulPartitionedCall�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
 dense_23/StatefulPartitionedCallStatefulPartitionedCalldense_23_inputdense_23_602610dense_23_602612*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_23_layer_call_and_return_conditional_losses_602487f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    ~
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602610*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602610*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602612*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602612*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_23/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_23/StatefulPartitionedCall-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2D
 dense_23/StatefulPartitionedCall dense_23/StatefulPartitionedCall2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:W S
'
_output_shapes
:���������
(
_user_specified_namedense_23_input
��
�#
"__inference__traced_restore_607720
file_prefix#
assignvariableop_variable: '
assignvariableop_1_variable_1: '
assignvariableop_2_variable_2: &
assignvariableop_3_adam_iter:	 (
assignvariableop_4_adam_beta_1: (
assignvariableop_5_adam_beta_2: '
assignvariableop_6_adam_decay: /
%assignvariableop_7_adam_learning_rate: 5
"assignvariableop_8_dense_16_kernel:	�/
 assignvariableop_9_dense_16_bias:	�7
#assignvariableop_10_dense_17_kernel:
��0
!assignvariableop_11_dense_17_bias:	�7
#assignvariableop_12_dense_18_kernel:
��0
!assignvariableop_13_dense_18_bias:	�7
#assignvariableop_14_dense_19_kernel:
��0
!assignvariableop_15_dense_19_bias:	�7
#assignvariableop_16_dense_20_kernel:
��0
!assignvariableop_17_dense_20_bias:	�7
#assignvariableop_18_dense_21_kernel:
��0
!assignvariableop_19_dense_21_bias:	�6
#assignvariableop_20_dense_22_kernel:	�/
!assignvariableop_21_dense_22_bias:5
#assignvariableop_22_dense_23_kernel:/
!assignvariableop_23_dense_23_bias:#
assignvariableop_24_total: #
assignvariableop_25_count: =
*assignvariableop_26_adam_dense_16_kernel_m:	�7
(assignvariableop_27_adam_dense_16_bias_m:	�>
*assignvariableop_28_adam_dense_17_kernel_m:
��7
(assignvariableop_29_adam_dense_17_bias_m:	�>
*assignvariableop_30_adam_dense_18_kernel_m:
��7
(assignvariableop_31_adam_dense_18_bias_m:	�>
*assignvariableop_32_adam_dense_19_kernel_m:
��7
(assignvariableop_33_adam_dense_19_bias_m:	�>
*assignvariableop_34_adam_dense_20_kernel_m:
��7
(assignvariableop_35_adam_dense_20_bias_m:	�>
*assignvariableop_36_adam_dense_21_kernel_m:
��7
(assignvariableop_37_adam_dense_21_bias_m:	�=
*assignvariableop_38_adam_dense_22_kernel_m:	�6
(assignvariableop_39_adam_dense_22_bias_m:<
*assignvariableop_40_adam_dense_23_kernel_m:6
(assignvariableop_41_adam_dense_23_bias_m:=
*assignvariableop_42_adam_dense_16_kernel_v:	�7
(assignvariableop_43_adam_dense_16_bias_v:	�>
*assignvariableop_44_adam_dense_17_kernel_v:
��7
(assignvariableop_45_adam_dense_17_bias_v:	�>
*assignvariableop_46_adam_dense_18_kernel_v:
��7
(assignvariableop_47_adam_dense_18_bias_v:	�>
*assignvariableop_48_adam_dense_19_kernel_v:
��7
(assignvariableop_49_adam_dense_19_bias_v:	�>
*assignvariableop_50_adam_dense_20_kernel_v:
��7
(assignvariableop_51_adam_dense_20_bias_v:	�>
*assignvariableop_52_adam_dense_21_kernel_v:
��7
(assignvariableop_53_adam_dense_21_bias_v:	�=
*assignvariableop_54_adam_dense_22_kernel_v:	�6
(assignvariableop_55_adam_dense_22_bias_v:<
*assignvariableop_56_adam_dense_23_kernel_v:6
(assignvariableop_57_adam_dense_23_bias_v:
identity_59��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_50�AssignVariableOp_51�AssignVariableOp_52�AssignVariableOp_53�AssignVariableOp_54�AssignVariableOp_55�AssignVariableOp_56�AssignVariableOp_57�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:;*
dtype0*�
value�B�;Ba1/.ATTRIBUTES/VARIABLE_VALUEBa2/.ATTRIBUTES/VARIABLE_VALUEBa3/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:;*
dtype0*�
value�B~;B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*I
dtypes?
=2;	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOpassignvariableop_variableIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOpassignvariableop_1_variable_1Identity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpassignvariableop_2_variable_2Identity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_iterIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpassignvariableop_4_adam_beta_1Identity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpassignvariableop_5_adam_beta_2Identity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_adam_decayIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOp%assignvariableop_7_adam_learning_rateIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp"assignvariableop_8_dense_16_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp assignvariableop_9_dense_16_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp#assignvariableop_10_dense_17_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOp!assignvariableop_11_dense_17_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp#assignvariableop_12_dense_18_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp!assignvariableop_13_dense_18_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp#assignvariableop_14_dense_19_kernelIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp!assignvariableop_15_dense_19_biasIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOp#assignvariableop_16_dense_20_kernelIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp!assignvariableop_17_dense_20_biasIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp#assignvariableop_18_dense_21_kernelIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp!assignvariableop_19_dense_21_biasIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOp#assignvariableop_20_dense_22_kernelIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOp!assignvariableop_21_dense_22_biasIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOp#assignvariableop_22_dense_23_kernelIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp!assignvariableop_23_dense_23_biasIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOpassignvariableop_24_totalIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpassignvariableop_25_countIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOp*assignvariableop_26_adam_dense_16_kernel_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOp(assignvariableop_27_adam_dense_16_bias_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOp*assignvariableop_28_adam_dense_17_kernel_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOp(assignvariableop_29_adam_dense_17_bias_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOp*assignvariableop_30_adam_dense_18_kernel_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOp(assignvariableop_31_adam_dense_18_bias_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOp*assignvariableop_32_adam_dense_19_kernel_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOp(assignvariableop_33_adam_dense_19_bias_mIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOp*assignvariableop_34_adam_dense_20_kernel_mIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOp(assignvariableop_35_adam_dense_20_bias_mIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOp*assignvariableop_36_adam_dense_21_kernel_mIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOp(assignvariableop_37_adam_dense_21_bias_mIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOp*assignvariableop_38_adam_dense_22_kernel_mIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOp(assignvariableop_39_adam_dense_22_bias_mIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOp*assignvariableop_40_adam_dense_23_kernel_mIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_41AssignVariableOp(assignvariableop_41_adam_dense_23_bias_mIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOp*assignvariableop_42_adam_dense_16_kernel_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_43AssignVariableOp(assignvariableop_43_adam_dense_16_bias_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_44AssignVariableOp*assignvariableop_44_adam_dense_17_kernel_vIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_45AssignVariableOp(assignvariableop_45_adam_dense_17_bias_vIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_46AssignVariableOp*assignvariableop_46_adam_dense_18_kernel_vIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_47AssignVariableOp(assignvariableop_47_adam_dense_18_bias_vIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_48AssignVariableOp*assignvariableop_48_adam_dense_19_kernel_vIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_49AssignVariableOp(assignvariableop_49_adam_dense_19_bias_vIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_50AssignVariableOp*assignvariableop_50_adam_dense_20_kernel_vIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_51AssignVariableOp(assignvariableop_51_adam_dense_20_bias_vIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_52AssignVariableOp*assignvariableop_52_adam_dense_21_kernel_vIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_53AssignVariableOp(assignvariableop_53_adam_dense_21_bias_vIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_54AssignVariableOp*assignvariableop_54_adam_dense_22_kernel_vIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_55AssignVariableOp(assignvariableop_55_adam_dense_22_bias_vIdentity_55:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_56AssignVariableOp*assignvariableop_56_adam_dense_23_kernel_vIdentity_56:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_57AssignVariableOp(assignvariableop_57_adam_dense_23_bias_vIdentity_57:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �

Identity_58Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_59IdentityIdentity_58:output:0^NoOp_1*
T0*
_output_shapes
: �

NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_59Identity_59:output:0*�
_input_shapesx
v: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
-__inference_sequential_4_layer_call_fn_601942
dense_16_input
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_16_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
'
_output_shapes
:���������
(
_user_specified_namedense_16_input
�
�
)__inference_dense_20_layer_call_fn_606739

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_20_layer_call_and_return_conditional_losses_601182p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�0
�
D__inference_dense_18_layer_call_and_return_conditional_losses_606620

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
)__inference_dense_17_layer_call_fn_606499

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_17_layer_call_and_return_conditional_losses_601041p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�0
�
D__inference_dense_22_layer_call_and_return_conditional_losses_601276

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
SeluSeluBiasAdd:output:0*
T0*'
_output_shapes
:���������f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentitySelu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�0
�
D__inference_dense_18_layer_call_and_return_conditional_losses_601088

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_0_606960J
7dense_16_kernel_regularizer_abs_readvariableop_resource:	�
identity��.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOpf
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_16_kernel_regularizer_abs_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_16_kernel_regularizer_abs_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_16/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp
�
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878

inputs"
dense_16_601632:	�
dense_16_601634:	�#
dense_17_601637:
��
dense_17_601639:	�#
dense_18_601642:
��
dense_18_601644:	�#
dense_19_601647:
��
dense_19_601649:	�#
dense_20_601652:
��
dense_20_601654:	�#
dense_21_601657:
��
dense_21_601659:	�"
dense_22_601662:	�
dense_22_601664:
identity�� dense_16/StatefulPartitionedCall�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp� dense_17/StatefulPartitionedCall�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp� dense_18/StatefulPartitionedCall�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp� dense_19/StatefulPartitionedCall�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp� dense_20/StatefulPartitionedCall�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp� dense_21/StatefulPartitionedCall�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp� dense_22/StatefulPartitionedCall�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
 dense_16/StatefulPartitionedCallStatefulPartitionedCallinputsdense_16_601632dense_16_601634*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_16_layer_call_and_return_conditional_losses_600994�
 dense_17/StatefulPartitionedCallStatefulPartitionedCall)dense_16/StatefulPartitionedCall:output:0dense_17_601637dense_17_601639*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_17_layer_call_and_return_conditional_losses_601041�
 dense_18/StatefulPartitionedCallStatefulPartitionedCall)dense_17/StatefulPartitionedCall:output:0dense_18_601642dense_18_601644*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_18_layer_call_and_return_conditional_losses_601088�
 dense_19/StatefulPartitionedCallStatefulPartitionedCall)dense_18/StatefulPartitionedCall:output:0dense_19_601647dense_19_601649*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_19_layer_call_and_return_conditional_losses_601135�
 dense_20/StatefulPartitionedCallStatefulPartitionedCall)dense_19/StatefulPartitionedCall:output:0dense_20_601652dense_20_601654*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_20_layer_call_and_return_conditional_losses_601182�
 dense_21/StatefulPartitionedCallStatefulPartitionedCall)dense_20/StatefulPartitionedCall:output:0dense_21_601657dense_21_601659*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_21_layer_call_and_return_conditional_losses_601229�
 dense_22/StatefulPartitionedCallStatefulPartitionedCall)dense_21/StatefulPartitionedCall:output:0dense_22_601662dense_22_601664*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_22_layer_call_and_return_conditional_losses_601276f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_601632*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_601632*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_16_601634*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_16_601634*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601637* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601637* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_17_601639*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_17_601639*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601642* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601642* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_18_601644*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_18_601644*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601647* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601647* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_19_601649*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_19_601649*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601652* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601652* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_20_601654*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_20_601654*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601657* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601657* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    y
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_21_601659*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: |
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_21_601659*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601662*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601662*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_22_601664*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_22_601664*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_22/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_16/StatefulPartitionedCall-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp!^dense_17/StatefulPartitionedCall-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp!^dense_18/StatefulPartitionedCall-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp!^dense_19/StatefulPartitionedCall-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp!^dense_20/StatefulPartitionedCall-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp!^dense_21/StatefulPartitionedCall-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp!^dense_22/StatefulPartitionedCall-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2D
 dense_16/StatefulPartitionedCall dense_16/StatefulPartitionedCall2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2D
 dense_17/StatefulPartitionedCall dense_17/StatefulPartitionedCall2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2D
 dense_18/StatefulPartitionedCall dense_18/StatefulPartitionedCall2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2D
 dense_19/StatefulPartitionedCall dense_19/StatefulPartitionedCall2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2D
 dense_20/StatefulPartitionedCall dense_20/StatefulPartitionedCall2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2D
 dense_21/StatefulPartitionedCall dense_21/StatefulPartitionedCall2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2D
 dense_22/StatefulPartitionedCall dense_22/StatefulPartitionedCall2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�0
�
D__inference_dense_16_layer_call_and_return_conditional_losses_600994

inputs1
matmul_readvariableop_resource:	�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_6_607080K
7dense_19_kernel_regularizer_abs_readvariableop_resource:
��
identity��.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOpf
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_19_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_19_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_19/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp
��
�
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_604279
input_1&
sequential_4_603930:	�"
sequential_4_603932:	�'
sequential_4_603934:
��"
sequential_4_603936:	�'
sequential_4_603938:
��"
sequential_4_603940:	�'
sequential_4_603942:
��"
sequential_4_603944:	�'
sequential_4_603946:
��"
sequential_4_603948:	�'
sequential_4_603950:
��"
sequential_4_603952:	�&
sequential_4_603954:	�!
sequential_4_603956:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: %
sequential_5_603981:!
sequential_5_603983:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�&sequential_4/StatefulPartitionedCall_1�$sequential_5/StatefulPartitionedCall�&sequential_5/StatefulPartitionedCall_1�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallinput_1sequential_4_603930sequential_4_603932sequential_4_603934sequential_4_603936sequential_4_603938sequential_4_603940sequential_4_603942sequential_4_603944sequential_4_603946sequential_4_603948sequential_4_603950sequential_4_603952sequential_4_603954sequential_4_603956*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
$sequential_5/StatefulPartitionedCallStatefulPartitionedCallstack:output:0sequential_5_603981sequential_5_603983*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591�
&sequential_5/StatefulPartitionedCall_1StatefulPartitionedCall-sequential_4/StatefulPartitionedCall:output:0sequential_5_603981sequential_5_603983*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591v
subSubinput_1/sequential_5/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSliceinput_1strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSliceinput_1strided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSliceinput_1strided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������
sub_2Sub-sequential_5/StatefulPartitionedCall:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCallstack_1:output:0sequential_4_603930sequential_4_603932sequential_4_603934sequential_4_603936sequential_4_603938sequential_4_603940sequential_4_603942sequential_4_603944sequential_4_603946sequential_4_603948sequential_4_603950sequential_4_603952sequential_4_603954sequential_4_603956*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878
sub_3Substack:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603930*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603930*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603932*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603932*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603934* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603934* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603936*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603936*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603938* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603938* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603940*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603940*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603942* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603942* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603944*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603944*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603946* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603946* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603948*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603948*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603950* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603950* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603952*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603952*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603954*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603954*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_603956*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_603956*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603981*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603981*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_603983*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_603983*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: |
IdentityIdentity-sequential_5/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1%^sequential_5/StatefulPartitionedCall'^sequential_5/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_12L
$sequential_5/StatefulPartitionedCall$sequential_5/StatefulPartitionedCall2P
&sequential_5/StatefulPartitionedCall_1&sequential_5/StatefulPartitionedCall_1:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�-
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591

inputs!
dense_23_602555:
dense_23_602557:
identity�� dense_23/StatefulPartitionedCall�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
 dense_23/StatefulPartitionedCallStatefulPartitionedCallinputsdense_23_602555dense_23_602557*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_23_layer_call_and_return_conditional_losses_602487f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    ~
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602555*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602555*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602557*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602557*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_23/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_23/StatefulPartitionedCall-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2D
 dense_23/StatefulPartitionedCall dense_23/StatefulPartitionedCall2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_10_607160K
7dense_21_kernel_regularizer_abs_readvariableop_resource:
��
identity��.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOpf
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp7dense_21_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_21_kernel_regularizer_abs_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: c
IdentityIdentity%dense_21/kernel/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp
��
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_606252

inputs:
'dense_16_matmul_readvariableop_resource:	�7
(dense_16_biasadd_readvariableop_resource:	�;
'dense_17_matmul_readvariableop_resource:
��7
(dense_17_biasadd_readvariableop_resource:	�;
'dense_18_matmul_readvariableop_resource:
��7
(dense_18_biasadd_readvariableop_resource:	�;
'dense_19_matmul_readvariableop_resource:
��7
(dense_19_biasadd_readvariableop_resource:	�;
'dense_20_matmul_readvariableop_resource:
��7
(dense_20_biasadd_readvariableop_resource:	�;
'dense_21_matmul_readvariableop_resource:
��7
(dense_21_biasadd_readvariableop_resource:	�:
'dense_22_matmul_readvariableop_resource:	�6
(dense_22_biasadd_readvariableop_resource:
identity��dense_16/BiasAdd/ReadVariableOp�dense_16/MatMul/ReadVariableOp�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�dense_17/BiasAdd/ReadVariableOp�dense_17/MatMul/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�dense_18/BiasAdd/ReadVariableOp�dense_18/MatMul/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�dense_19/BiasAdd/ReadVariableOp�dense_19/MatMul/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�dense_20/BiasAdd/ReadVariableOp�dense_20/MatMul/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�dense_21/BiasAdd/ReadVariableOp�dense_21/MatMul/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�dense_22/BiasAdd/ReadVariableOp�dense_22/MatMul/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
dense_16/MatMul/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0|
dense_16/MatMulMatMulinputs&dense_16/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_16/BiasAdd/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/BiasAddBiasAdddense_16/MatMul:product:0'dense_16/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_16/SeluSeludense_16/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_17/MatMul/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/MatMulMatMuldense_16/Selu:activations:0&dense_17/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_17/BiasAdd/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/BiasAddBiasAdddense_17/MatMul:product:0'dense_17/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_17/SeluSeludense_17/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_18/MatMul/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/MatMulMatMuldense_17/Selu:activations:0&dense_18/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_18/BiasAdd/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/BiasAddBiasAdddense_18/MatMul:product:0'dense_18/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_18/SeluSeludense_18/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_19/MatMul/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/MatMulMatMuldense_18/Selu:activations:0&dense_19/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_19/BiasAdd/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/BiasAddBiasAdddense_19/MatMul:product:0'dense_19/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_19/SeluSeludense_19/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_20/MatMul/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/MatMulMatMuldense_19/Selu:activations:0&dense_20/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_20/BiasAdd/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/BiasAddBiasAdddense_20/MatMul:product:0'dense_20/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_20/SeluSeludense_20/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_21/MatMul/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/MatMulMatMuldense_20/Selu:activations:0&dense_21/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_21/BiasAdd/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/BiasAddBiasAdddense_21/MatMul:product:0'dense_21/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_21/SeluSeludense_21/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_22/MatMul/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/MatMulMatMuldense_21/Selu:activations:0&dense_22/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
dense_22/BiasAdd/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_22/BiasAddBiasAdddense_22/MatMul:product:0'dense_22/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_22/SeluSeludense_22/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: j
IdentityIdentitydense_22/Selu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp ^dense_16/BiasAdd/ReadVariableOp^dense_16/MatMul/ReadVariableOp-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp ^dense_17/BiasAdd/ReadVariableOp^dense_17/MatMul/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp ^dense_18/BiasAdd/ReadVariableOp^dense_18/MatMul/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp ^dense_19/BiasAdd/ReadVariableOp^dense_19/MatMul/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp ^dense_20/BiasAdd/ReadVariableOp^dense_20/MatMul/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp ^dense_21/BiasAdd/ReadVariableOp^dense_21/MatMul/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp ^dense_22/BiasAdd/ReadVariableOp^dense_22/MatMul/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2B
dense_16/BiasAdd/ReadVariableOpdense_16/BiasAdd/ReadVariableOp2@
dense_16/MatMul/ReadVariableOpdense_16/MatMul/ReadVariableOp2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2B
dense_17/BiasAdd/ReadVariableOpdense_17/BiasAdd/ReadVariableOp2@
dense_17/MatMul/ReadVariableOpdense_17/MatMul/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2B
dense_18/BiasAdd/ReadVariableOpdense_18/BiasAdd/ReadVariableOp2@
dense_18/MatMul/ReadVariableOpdense_18/MatMul/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2B
dense_19/BiasAdd/ReadVariableOpdense_19/BiasAdd/ReadVariableOp2@
dense_19/MatMul/ReadVariableOpdense_19/MatMul/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2B
dense_20/BiasAdd/ReadVariableOpdense_20/BiasAdd/ReadVariableOp2@
dense_20/MatMul/ReadVariableOpdense_20/MatMul/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2B
dense_21/BiasAdd/ReadVariableOpdense_21/BiasAdd/ReadVariableOp2@
dense_21/MatMul/ReadVariableOpdense_21/MatMul/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2B
dense_22/BiasAdd/ReadVariableOpdense_22/BiasAdd/ReadVariableOp2@
dense_22/MatMul/ReadVariableOpdense_22/MatMul/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
$__inference_signature_wrapper_604570
input_1
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:

unknown_13: 

unknown_14: 

unknown_15: 

unknown_16:

unknown_17:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*5
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__wrapped_model_600946o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
'
_output_shapes
:���������
!
_user_specified_name	input_1
�0
�
D__inference_dense_19_layer_call_and_return_conditional_losses_601135

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_3_607020D
5dense_17_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOpd
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_17_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_17_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_17/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp
�
�
)__inference_dense_23_layer_call_fn_607259

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_23_layer_call_and_return_conditional_losses_602487o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
܍
�#
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605450
xG
4sequential_4_dense_16_matmul_readvariableop_resource:	�D
5sequential_4_dense_16_biasadd_readvariableop_resource:	�H
4sequential_4_dense_17_matmul_readvariableop_resource:
��D
5sequential_4_dense_17_biasadd_readvariableop_resource:	�H
4sequential_4_dense_18_matmul_readvariableop_resource:
��D
5sequential_4_dense_18_biasadd_readvariableop_resource:	�H
4sequential_4_dense_19_matmul_readvariableop_resource:
��D
5sequential_4_dense_19_biasadd_readvariableop_resource:	�H
4sequential_4_dense_20_matmul_readvariableop_resource:
��D
5sequential_4_dense_20_biasadd_readvariableop_resource:	�H
4sequential_4_dense_21_matmul_readvariableop_resource:
��D
5sequential_4_dense_21_biasadd_readvariableop_resource:	�G
4sequential_4_dense_22_matmul_readvariableop_resource:	�C
5sequential_4_dense_22_biasadd_readvariableop_resource:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: F
4sequential_5_dense_23_matmul_readvariableop_resource:C
5sequential_5_dense_23_biasadd_readvariableop_resource:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�,sequential_4/dense_16/BiasAdd/ReadVariableOp�.sequential_4/dense_16/BiasAdd_1/ReadVariableOp�+sequential_4/dense_16/MatMul/ReadVariableOp�-sequential_4/dense_16/MatMul_1/ReadVariableOp�,sequential_4/dense_17/BiasAdd/ReadVariableOp�.sequential_4/dense_17/BiasAdd_1/ReadVariableOp�+sequential_4/dense_17/MatMul/ReadVariableOp�-sequential_4/dense_17/MatMul_1/ReadVariableOp�,sequential_4/dense_18/BiasAdd/ReadVariableOp�.sequential_4/dense_18/BiasAdd_1/ReadVariableOp�+sequential_4/dense_18/MatMul/ReadVariableOp�-sequential_4/dense_18/MatMul_1/ReadVariableOp�,sequential_4/dense_19/BiasAdd/ReadVariableOp�.sequential_4/dense_19/BiasAdd_1/ReadVariableOp�+sequential_4/dense_19/MatMul/ReadVariableOp�-sequential_4/dense_19/MatMul_1/ReadVariableOp�,sequential_4/dense_20/BiasAdd/ReadVariableOp�.sequential_4/dense_20/BiasAdd_1/ReadVariableOp�+sequential_4/dense_20/MatMul/ReadVariableOp�-sequential_4/dense_20/MatMul_1/ReadVariableOp�,sequential_4/dense_21/BiasAdd/ReadVariableOp�.sequential_4/dense_21/BiasAdd_1/ReadVariableOp�+sequential_4/dense_21/MatMul/ReadVariableOp�-sequential_4/dense_21/MatMul_1/ReadVariableOp�,sequential_4/dense_22/BiasAdd/ReadVariableOp�.sequential_4/dense_22/BiasAdd_1/ReadVariableOp�+sequential_4/dense_22/MatMul/ReadVariableOp�-sequential_4/dense_22/MatMul_1/ReadVariableOp�,sequential_5/dense_23/BiasAdd/ReadVariableOp�.sequential_5/dense_23/BiasAdd_1/ReadVariableOp�+sequential_5/dense_23/MatMul/ReadVariableOp�-sequential_5/dense_23/MatMul_1/ReadVariableOp�
+sequential_4/dense_16/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_16/MatMulMatMulx3sequential_4/dense_16/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_16/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_16/BiasAddBiasAdd&sequential_4/dense_16/MatMul:product:04sequential_4/dense_16/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_16/SeluSelu&sequential_4/dense_16/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_17/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_17/MatMulMatMul(sequential_4/dense_16/Selu:activations:03sequential_4/dense_17/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_17/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_17/BiasAddBiasAdd&sequential_4/dense_17/MatMul:product:04sequential_4/dense_17/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_17/SeluSelu&sequential_4/dense_17/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_18/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_18/MatMulMatMul(sequential_4/dense_17/Selu:activations:03sequential_4/dense_18/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_18/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_18/BiasAddBiasAdd&sequential_4/dense_18/MatMul:product:04sequential_4/dense_18/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_18/SeluSelu&sequential_4/dense_18/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_19/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_19/MatMulMatMul(sequential_4/dense_18/Selu:activations:03sequential_4/dense_19/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_19/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_19/BiasAddBiasAdd&sequential_4/dense_19/MatMul:product:04sequential_4/dense_19/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_19/SeluSelu&sequential_4/dense_19/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_20/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_20/MatMulMatMul(sequential_4/dense_19/Selu:activations:03sequential_4/dense_20/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_20/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_20/BiasAddBiasAdd&sequential_4/dense_20/MatMul:product:04sequential_4/dense_20/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_20/SeluSelu&sequential_4/dense_20/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_21/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_21/MatMulMatMul(sequential_4/dense_20/Selu:activations:03sequential_4/dense_21/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
,sequential_4/dense_21/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_21/BiasAddBiasAdd&sequential_4/dense_21/MatMul:product:04sequential_4/dense_21/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������}
sequential_4/dense_21/SeluSelu&sequential_4/dense_21/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
+sequential_4/dense_22/MatMul/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_22/MatMulMatMul(sequential_4/dense_21/Selu:activations:03sequential_4/dense_22/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
,sequential_4/dense_22/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_4/dense_22/BiasAddBiasAdd&sequential_4/dense_22/MatMul:product:04sequential_4/dense_22/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������|
sequential_4/dense_22/SeluSelu&sequential_4/dense_22/BiasAdd:output:0*
T0*'
_output_shapes
:���������d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice(sequential_4/dense_22/Selu:activations:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
+sequential_5/dense_23/MatMul/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sequential_5/dense_23/MatMulMatMulstack:output:03sequential_5/dense_23/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
,sequential_5/dense_23/BiasAdd/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_5/dense_23/BiasAddBiasAdd&sequential_5/dense_23/MatMul:product:04sequential_5/dense_23/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
-sequential_5/dense_23/MatMul_1/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
sequential_5/dense_23/MatMul_1MatMul(sequential_4/dense_22/Selu:activations:05sequential_5/dense_23/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_5/dense_23/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_5/dense_23/BiasAdd_1BiasAdd(sequential_5/dense_23/MatMul_1:product:06sequential_5/dense_23/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������i
subSubx(sequential_5/dense_23/BiasAdd_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSlicexstrided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSlicexstrided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSlicexstrided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������x
sub_2Sub&sequential_5/dense_23/BiasAdd:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
-sequential_4/dense_16/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_16/MatMul_1MatMulstack_1:output:05sequential_4/dense_16/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_16/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_16/BiasAdd_1BiasAdd(sequential_4/dense_16/MatMul_1:product:06sequential_4/dense_16/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_16/Selu_1Selu(sequential_4/dense_16/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_17/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_17/MatMul_1MatMul*sequential_4/dense_16/Selu_1:activations:05sequential_4/dense_17/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_17/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_17/BiasAdd_1BiasAdd(sequential_4/dense_17/MatMul_1:product:06sequential_4/dense_17/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_17/Selu_1Selu(sequential_4/dense_17/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_18/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_18/MatMul_1MatMul*sequential_4/dense_17/Selu_1:activations:05sequential_4/dense_18/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_18/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_18/BiasAdd_1BiasAdd(sequential_4/dense_18/MatMul_1:product:06sequential_4/dense_18/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_18/Selu_1Selu(sequential_4/dense_18/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_19/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_19/MatMul_1MatMul*sequential_4/dense_18/Selu_1:activations:05sequential_4/dense_19/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_19/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_19/BiasAdd_1BiasAdd(sequential_4/dense_19/MatMul_1:product:06sequential_4/dense_19/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_19/Selu_1Selu(sequential_4/dense_19/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_20/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_20/MatMul_1MatMul*sequential_4/dense_19/Selu_1:activations:05sequential_4/dense_20/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_20/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_20/BiasAdd_1BiasAdd(sequential_4/dense_20/MatMul_1:product:06sequential_4/dense_20/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_20/Selu_1Selu(sequential_4/dense_20/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_21/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_4/dense_21/MatMul_1MatMul*sequential_4/dense_20/Selu_1:activations:05sequential_4/dense_21/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
.sequential_4/dense_21/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
sequential_4/dense_21/BiasAdd_1BiasAdd(sequential_4/dense_21/MatMul_1:product:06sequential_4/dense_21/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_4/dense_21/Selu_1Selu(sequential_4/dense_21/BiasAdd_1:output:0*
T0*(
_output_shapes
:�����������
-sequential_4/dense_22/MatMul_1/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_4/dense_22/MatMul_1MatMul*sequential_4/dense_21/Selu_1:activations:05sequential_4/dense_22/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_4/dense_22/BiasAdd_1/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_4/dense_22/BiasAdd_1BiasAdd(sequential_4/dense_22/MatMul_1:product:06sequential_4/dense_22/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
sequential_4/dense_22/Selu_1Selu(sequential_4/dense_22/BiasAdd_1:output:0*
T0*'
_output_shapes
:���������z
sub_3Substack:output:0*sequential_4/dense_22/Selu_1:activations:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_4_dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_4_dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOp4sequential_5_dense_23_matmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOp5sequential_5_dense_23_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: u
IdentityIdentity&sequential_5/dense_23/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp-^sequential_4/dense_16/BiasAdd/ReadVariableOp/^sequential_4/dense_16/BiasAdd_1/ReadVariableOp,^sequential_4/dense_16/MatMul/ReadVariableOp.^sequential_4/dense_16/MatMul_1/ReadVariableOp-^sequential_4/dense_17/BiasAdd/ReadVariableOp/^sequential_4/dense_17/BiasAdd_1/ReadVariableOp,^sequential_4/dense_17/MatMul/ReadVariableOp.^sequential_4/dense_17/MatMul_1/ReadVariableOp-^sequential_4/dense_18/BiasAdd/ReadVariableOp/^sequential_4/dense_18/BiasAdd_1/ReadVariableOp,^sequential_4/dense_18/MatMul/ReadVariableOp.^sequential_4/dense_18/MatMul_1/ReadVariableOp-^sequential_4/dense_19/BiasAdd/ReadVariableOp/^sequential_4/dense_19/BiasAdd_1/ReadVariableOp,^sequential_4/dense_19/MatMul/ReadVariableOp.^sequential_4/dense_19/MatMul_1/ReadVariableOp-^sequential_4/dense_20/BiasAdd/ReadVariableOp/^sequential_4/dense_20/BiasAdd_1/ReadVariableOp,^sequential_4/dense_20/MatMul/ReadVariableOp.^sequential_4/dense_20/MatMul_1/ReadVariableOp-^sequential_4/dense_21/BiasAdd/ReadVariableOp/^sequential_4/dense_21/BiasAdd_1/ReadVariableOp,^sequential_4/dense_21/MatMul/ReadVariableOp.^sequential_4/dense_21/MatMul_1/ReadVariableOp-^sequential_4/dense_22/BiasAdd/ReadVariableOp/^sequential_4/dense_22/BiasAdd_1/ReadVariableOp,^sequential_4/dense_22/MatMul/ReadVariableOp.^sequential_4/dense_22/MatMul_1/ReadVariableOp-^sequential_5/dense_23/BiasAdd/ReadVariableOp/^sequential_5/dense_23/BiasAdd_1/ReadVariableOp,^sequential_5/dense_23/MatMul/ReadVariableOp.^sequential_5/dense_23/MatMul_1/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2\
,sequential_4/dense_16/BiasAdd/ReadVariableOp,sequential_4/dense_16/BiasAdd/ReadVariableOp2`
.sequential_4/dense_16/BiasAdd_1/ReadVariableOp.sequential_4/dense_16/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_16/MatMul/ReadVariableOp+sequential_4/dense_16/MatMul/ReadVariableOp2^
-sequential_4/dense_16/MatMul_1/ReadVariableOp-sequential_4/dense_16/MatMul_1/ReadVariableOp2\
,sequential_4/dense_17/BiasAdd/ReadVariableOp,sequential_4/dense_17/BiasAdd/ReadVariableOp2`
.sequential_4/dense_17/BiasAdd_1/ReadVariableOp.sequential_4/dense_17/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_17/MatMul/ReadVariableOp+sequential_4/dense_17/MatMul/ReadVariableOp2^
-sequential_4/dense_17/MatMul_1/ReadVariableOp-sequential_4/dense_17/MatMul_1/ReadVariableOp2\
,sequential_4/dense_18/BiasAdd/ReadVariableOp,sequential_4/dense_18/BiasAdd/ReadVariableOp2`
.sequential_4/dense_18/BiasAdd_1/ReadVariableOp.sequential_4/dense_18/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_18/MatMul/ReadVariableOp+sequential_4/dense_18/MatMul/ReadVariableOp2^
-sequential_4/dense_18/MatMul_1/ReadVariableOp-sequential_4/dense_18/MatMul_1/ReadVariableOp2\
,sequential_4/dense_19/BiasAdd/ReadVariableOp,sequential_4/dense_19/BiasAdd/ReadVariableOp2`
.sequential_4/dense_19/BiasAdd_1/ReadVariableOp.sequential_4/dense_19/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_19/MatMul/ReadVariableOp+sequential_4/dense_19/MatMul/ReadVariableOp2^
-sequential_4/dense_19/MatMul_1/ReadVariableOp-sequential_4/dense_19/MatMul_1/ReadVariableOp2\
,sequential_4/dense_20/BiasAdd/ReadVariableOp,sequential_4/dense_20/BiasAdd/ReadVariableOp2`
.sequential_4/dense_20/BiasAdd_1/ReadVariableOp.sequential_4/dense_20/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_20/MatMul/ReadVariableOp+sequential_4/dense_20/MatMul/ReadVariableOp2^
-sequential_4/dense_20/MatMul_1/ReadVariableOp-sequential_4/dense_20/MatMul_1/ReadVariableOp2\
,sequential_4/dense_21/BiasAdd/ReadVariableOp,sequential_4/dense_21/BiasAdd/ReadVariableOp2`
.sequential_4/dense_21/BiasAdd_1/ReadVariableOp.sequential_4/dense_21/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_21/MatMul/ReadVariableOp+sequential_4/dense_21/MatMul/ReadVariableOp2^
-sequential_4/dense_21/MatMul_1/ReadVariableOp-sequential_4/dense_21/MatMul_1/ReadVariableOp2\
,sequential_4/dense_22/BiasAdd/ReadVariableOp,sequential_4/dense_22/BiasAdd/ReadVariableOp2`
.sequential_4/dense_22/BiasAdd_1/ReadVariableOp.sequential_4/dense_22/BiasAdd_1/ReadVariableOp2Z
+sequential_4/dense_22/MatMul/ReadVariableOp+sequential_4/dense_22/MatMul/ReadVariableOp2^
-sequential_4/dense_22/MatMul_1/ReadVariableOp-sequential_4/dense_22/MatMul_1/ReadVariableOp2\
,sequential_5/dense_23/BiasAdd/ReadVariableOp,sequential_5/dense_23/BiasAdd/ReadVariableOp2`
.sequential_5/dense_23/BiasAdd_1/ReadVariableOp.sequential_5/dense_23/BiasAdd_1/ReadVariableOp2Z
+sequential_5/dense_23/MatMul/ReadVariableOp+sequential_5/dense_23/MatMul/ReadVariableOp2^
-sequential_5/dense_23/MatMul_1/ReadVariableOp-sequential_5/dense_23/MatMul_1/ReadVariableOp:J F
'
_output_shapes
:���������

_user_specified_namex
�
�
)__inference_dense_19_layer_call_fn_606659

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_19_layer_call_and_return_conditional_losses_601135p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
-__inference_sequential_5_layer_call_fn_606300

inputs
unknown:
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602591o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�.
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_602685
dense_23_input!
dense_23_602649:
dense_23_602651:
identity�� dense_23/StatefulPartitionedCall�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
 dense_23/StatefulPartitionedCallStatefulPartitionedCalldense_23_inputdense_23_602649dense_23_602651*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_23_layer_call_and_return_conditional_losses_602487f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    ~
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602649*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602649*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602651*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602651*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_23/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_23/StatefulPartitionedCall-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2D
 dense_23/StatefulPartitionedCall dense_23/StatefulPartitionedCall2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:W S
'
_output_shapes
:���������
(
_user_specified_namedense_23_input
�
�
)__inference_dense_18_layer_call_fn_606579

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_18_layer_call_and_return_conditional_losses_601088p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�/
�
D__inference_dense_23_layer_call_and_return_conditional_losses_607299

inputs0
matmul_readvariableop_resource:-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: _
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�0
�
D__inference_dense_21_layer_call_and_return_conditional_losses_606860

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�k
�
__inference__traced_save_607536
file_prefix'
#savev2_variable_read_readvariableop)
%savev2_variable_1_read_readvariableop)
%savev2_variable_2_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop.
*savev2_dense_16_kernel_read_readvariableop,
(savev2_dense_16_bias_read_readvariableop.
*savev2_dense_17_kernel_read_readvariableop,
(savev2_dense_17_bias_read_readvariableop.
*savev2_dense_18_kernel_read_readvariableop,
(savev2_dense_18_bias_read_readvariableop.
*savev2_dense_19_kernel_read_readvariableop,
(savev2_dense_19_bias_read_readvariableop.
*savev2_dense_20_kernel_read_readvariableop,
(savev2_dense_20_bias_read_readvariableop.
*savev2_dense_21_kernel_read_readvariableop,
(savev2_dense_21_bias_read_readvariableop.
*savev2_dense_22_kernel_read_readvariableop,
(savev2_dense_22_bias_read_readvariableop.
*savev2_dense_23_kernel_read_readvariableop,
(savev2_dense_23_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop5
1savev2_adam_dense_16_kernel_m_read_readvariableop3
/savev2_adam_dense_16_bias_m_read_readvariableop5
1savev2_adam_dense_17_kernel_m_read_readvariableop3
/savev2_adam_dense_17_bias_m_read_readvariableop5
1savev2_adam_dense_18_kernel_m_read_readvariableop3
/savev2_adam_dense_18_bias_m_read_readvariableop5
1savev2_adam_dense_19_kernel_m_read_readvariableop3
/savev2_adam_dense_19_bias_m_read_readvariableop5
1savev2_adam_dense_20_kernel_m_read_readvariableop3
/savev2_adam_dense_20_bias_m_read_readvariableop5
1savev2_adam_dense_21_kernel_m_read_readvariableop3
/savev2_adam_dense_21_bias_m_read_readvariableop5
1savev2_adam_dense_22_kernel_m_read_readvariableop3
/savev2_adam_dense_22_bias_m_read_readvariableop5
1savev2_adam_dense_23_kernel_m_read_readvariableop3
/savev2_adam_dense_23_bias_m_read_readvariableop5
1savev2_adam_dense_16_kernel_v_read_readvariableop3
/savev2_adam_dense_16_bias_v_read_readvariableop5
1savev2_adam_dense_17_kernel_v_read_readvariableop3
/savev2_adam_dense_17_bias_v_read_readvariableop5
1savev2_adam_dense_18_kernel_v_read_readvariableop3
/savev2_adam_dense_18_bias_v_read_readvariableop5
1savev2_adam_dense_19_kernel_v_read_readvariableop3
/savev2_adam_dense_19_bias_v_read_readvariableop5
1savev2_adam_dense_20_kernel_v_read_readvariableop3
/savev2_adam_dense_20_bias_v_read_readvariableop5
1savev2_adam_dense_21_kernel_v_read_readvariableop3
/savev2_adam_dense_21_bias_v_read_readvariableop5
1savev2_adam_dense_22_kernel_v_read_readvariableop3
/savev2_adam_dense_22_bias_v_read_readvariableop5
1savev2_adam_dense_23_kernel_v_read_readvariableop3
/savev2_adam_dense_23_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:;*
dtype0*�
value�B�;Ba1/.ATTRIBUTES/VARIABLE_VALUEBa2/.ATTRIBUTES/VARIABLE_VALUEBa3/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB&variables/6/.ATTRIBUTES/VARIABLE_VALUEB&variables/7/.ATTRIBUTES/VARIABLE_VALUEB&variables/8/.ATTRIBUTES/VARIABLE_VALUEB&variables/9/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/12/.ATTRIBUTES/VARIABLE_VALUEB'variables/13/.ATTRIBUTES/VARIABLE_VALUEB'variables/14/.ATTRIBUTES/VARIABLE_VALUEB'variables/15/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBBvariables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBBvariables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBCvariables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:;*
dtype0*�
value�B~;B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0#savev2_variable_read_readvariableop%savev2_variable_1_read_readvariableop%savev2_variable_2_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop*savev2_dense_16_kernel_read_readvariableop(savev2_dense_16_bias_read_readvariableop*savev2_dense_17_kernel_read_readvariableop(savev2_dense_17_bias_read_readvariableop*savev2_dense_18_kernel_read_readvariableop(savev2_dense_18_bias_read_readvariableop*savev2_dense_19_kernel_read_readvariableop(savev2_dense_19_bias_read_readvariableop*savev2_dense_20_kernel_read_readvariableop(savev2_dense_20_bias_read_readvariableop*savev2_dense_21_kernel_read_readvariableop(savev2_dense_21_bias_read_readvariableop*savev2_dense_22_kernel_read_readvariableop(savev2_dense_22_bias_read_readvariableop*savev2_dense_23_kernel_read_readvariableop(savev2_dense_23_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop1savev2_adam_dense_16_kernel_m_read_readvariableop/savev2_adam_dense_16_bias_m_read_readvariableop1savev2_adam_dense_17_kernel_m_read_readvariableop/savev2_adam_dense_17_bias_m_read_readvariableop1savev2_adam_dense_18_kernel_m_read_readvariableop/savev2_adam_dense_18_bias_m_read_readvariableop1savev2_adam_dense_19_kernel_m_read_readvariableop/savev2_adam_dense_19_bias_m_read_readvariableop1savev2_adam_dense_20_kernel_m_read_readvariableop/savev2_adam_dense_20_bias_m_read_readvariableop1savev2_adam_dense_21_kernel_m_read_readvariableop/savev2_adam_dense_21_bias_m_read_readvariableop1savev2_adam_dense_22_kernel_m_read_readvariableop/savev2_adam_dense_22_bias_m_read_readvariableop1savev2_adam_dense_23_kernel_m_read_readvariableop/savev2_adam_dense_23_bias_m_read_readvariableop1savev2_adam_dense_16_kernel_v_read_readvariableop/savev2_adam_dense_16_bias_v_read_readvariableop1savev2_adam_dense_17_kernel_v_read_readvariableop/savev2_adam_dense_17_bias_v_read_readvariableop1savev2_adam_dense_18_kernel_v_read_readvariableop/savev2_adam_dense_18_bias_v_read_readvariableop1savev2_adam_dense_19_kernel_v_read_readvariableop/savev2_adam_dense_19_bias_v_read_readvariableop1savev2_adam_dense_20_kernel_v_read_readvariableop/savev2_adam_dense_20_bias_v_read_readvariableop1savev2_adam_dense_21_kernel_v_read_readvariableop/savev2_adam_dense_21_bias_v_read_readvariableop1savev2_adam_dense_22_kernel_v_read_readvariableop/savev2_adam_dense_22_bias_v_read_readvariableop1savev2_adam_dense_23_kernel_v_read_readvariableop/savev2_adam_dense_23_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *I
dtypes?
=2;	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*�
_input_shapes�
�: : : : : : : : : :	�:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:::: : :	�:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�::::	�:�:
��:�:
��:�:
��:�:
��:�:
��:�:	�:::: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :%	!

_output_shapes
:	�:!


_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%!

_output_shapes
:	�: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	�:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:! 

_output_shapes	
:�:&!"
 
_output_shapes
:
��:!"

_output_shapes	
:�:&#"
 
_output_shapes
:
��:!$

_output_shapes	
:�:&%"
 
_output_shapes
:
��:!&

_output_shapes	
:�:%'!

_output_shapes
:	�: (

_output_shapes
::$) 

_output_shapes

:: *

_output_shapes
::%+!

_output_shapes
:	�:!,

_output_shapes	
:�:&-"
 
_output_shapes
:
��:!.

_output_shapes	
:�:&/"
 
_output_shapes
:
��:!0

_output_shapes	
:�:&1"
 
_output_shapes
:
��:!2

_output_shapes	
:�:&3"
 
_output_shapes
:
��:!4

_output_shapes	
:�:&5"
 
_output_shapes
:
��:!6

_output_shapes	
:�:%7!

_output_shapes
:	�: 8

_output_shapes
::$9 

_output_shapes

:: :

_output_shapes
::;

_output_shapes
: 
�-
�
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524

inputs!
dense_23_602488:
dense_23_602490:
identity�� dense_23/StatefulPartitionedCall�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�
 dense_23/StatefulPartitionedCallStatefulPartitionedCallinputsdense_23_602488dense_23_602490*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_23_layer_call_and_return_conditional_losses_602487f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    ~
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602488*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602488*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    x
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpdense_23_602490*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: {
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_23_602490*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: x
IdentityIdentity)dense_23/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp!^dense_23/StatefulPartitionedCall-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 2D
 dense_23/StatefulPartitionedCall dense_23/StatefulPartitionedCall2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
H__inference_sequential_4_layer_call_and_return_conditional_losses_605989

inputs:
'dense_16_matmul_readvariableop_resource:	�7
(dense_16_biasadd_readvariableop_resource:	�;
'dense_17_matmul_readvariableop_resource:
��7
(dense_17_biasadd_readvariableop_resource:	�;
'dense_18_matmul_readvariableop_resource:
��7
(dense_18_biasadd_readvariableop_resource:	�;
'dense_19_matmul_readvariableop_resource:
��7
(dense_19_biasadd_readvariableop_resource:	�;
'dense_20_matmul_readvariableop_resource:
��7
(dense_20_biasadd_readvariableop_resource:	�;
'dense_21_matmul_readvariableop_resource:
��7
(dense_21_biasadd_readvariableop_resource:	�:
'dense_22_matmul_readvariableop_resource:	�6
(dense_22_biasadd_readvariableop_resource:
identity��dense_16/BiasAdd/ReadVariableOp�dense_16/MatMul/ReadVariableOp�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�dense_17/BiasAdd/ReadVariableOp�dense_17/MatMul/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�dense_18/BiasAdd/ReadVariableOp�dense_18/MatMul/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�dense_19/BiasAdd/ReadVariableOp�dense_19/MatMul/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�dense_20/BiasAdd/ReadVariableOp�dense_20/MatMul/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�dense_21/BiasAdd/ReadVariableOp�dense_21/MatMul/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�dense_22/BiasAdd/ReadVariableOp�dense_22/MatMul/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�
dense_16/MatMul/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0|
dense_16/MatMulMatMulinputs&dense_16/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_16/BiasAdd/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/BiasAddBiasAdddense_16/MatMul:product:0'dense_16/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_16/SeluSeludense_16/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_17/MatMul/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/MatMulMatMuldense_16/Selu:activations:0&dense_17/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_17/BiasAdd/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/BiasAddBiasAdddense_17/MatMul:product:0'dense_17/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_17/SeluSeludense_17/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_18/MatMul/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/MatMulMatMuldense_17/Selu:activations:0&dense_18/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_18/BiasAdd/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/BiasAddBiasAdddense_18/MatMul:product:0'dense_18/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_18/SeluSeludense_18/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_19/MatMul/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/MatMulMatMuldense_18/Selu:activations:0&dense_19/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_19/BiasAdd/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/BiasAddBiasAdddense_19/MatMul:product:0'dense_19/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_19/SeluSeludense_19/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_20/MatMul/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/MatMulMatMuldense_19/Selu:activations:0&dense_20/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_20/BiasAdd/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/BiasAddBiasAdddense_20/MatMul:product:0'dense_20/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_20/SeluSeludense_20/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_21/MatMul/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/MatMulMatMuldense_20/Selu:activations:0&dense_21/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
dense_21/BiasAdd/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/BiasAddBiasAdddense_21/MatMul:product:0'dense_21/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������c
dense_21/SeluSeludense_21/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
dense_22/MatMul/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/MatMulMatMuldense_21/Selu:activations:0&dense_22/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
dense_22/BiasAdd/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
dense_22/BiasAddBiasAdddense_22/MatMul:product:0'dense_22/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������b
dense_22/SeluSeludense_22/BiasAdd:output:0*
T0*'
_output_shapes
:���������f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_16_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_16_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_17_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_17_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_18_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_18_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_19_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_19_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_20_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_20_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_21_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_21_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'dense_22_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOp(dense_22_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: j
IdentityIdentitydense_22/Selu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp ^dense_16/BiasAdd/ReadVariableOp^dense_16/MatMul/ReadVariableOp-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp ^dense_17/BiasAdd/ReadVariableOp^dense_17/MatMul/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp ^dense_18/BiasAdd/ReadVariableOp^dense_18/MatMul/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp ^dense_19/BiasAdd/ReadVariableOp^dense_19/MatMul/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp ^dense_20/BiasAdd/ReadVariableOp^dense_20/MatMul/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp ^dense_21/BiasAdd/ReadVariableOp^dense_21/MatMul/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp ^dense_22/BiasAdd/ReadVariableOp^dense_22/MatMul/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 2B
dense_16/BiasAdd/ReadVariableOpdense_16/BiasAdd/ReadVariableOp2@
dense_16/MatMul/ReadVariableOpdense_16/MatMul/ReadVariableOp2\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2B
dense_17/BiasAdd/ReadVariableOpdense_17/BiasAdd/ReadVariableOp2@
dense_17/MatMul/ReadVariableOpdense_17/MatMul/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2B
dense_18/BiasAdd/ReadVariableOpdense_18/BiasAdd/ReadVariableOp2@
dense_18/MatMul/ReadVariableOpdense_18/MatMul/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2B
dense_19/BiasAdd/ReadVariableOpdense_19/BiasAdd/ReadVariableOp2@
dense_19/MatMul/ReadVariableOpdense_19/MatMul/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2B
dense_20/BiasAdd/ReadVariableOpdense_20/BiasAdd/ReadVariableOp2@
dense_20/MatMul/ReadVariableOpdense_20/MatMul/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2B
dense_21/BiasAdd/ReadVariableOpdense_21/BiasAdd/ReadVariableOp2@
dense_21/MatMul/ReadVariableOpdense_21/MatMul/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2B
dense_22/BiasAdd/ReadVariableOpdense_22/BiasAdd/ReadVariableOp2@
dense_22/MatMul/ReadVariableOpdense_22/MatMul/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_5_607060D
5dense_18_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOpd
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_18_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_18_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_18/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp
�0
�
D__inference_dense_22_layer_call_and_return_conditional_losses_606940

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
SeluSeluBiasAdd:output:0*
T0*'
_output_shapes
:���������f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentitySelu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
݁
�
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603041
x&
sequential_4_602692:	�"
sequential_4_602694:	�'
sequential_4_602696:
��"
sequential_4_602698:	�'
sequential_4_602700:
��"
sequential_4_602702:	�'
sequential_4_602704:
��"
sequential_4_602706:	�'
sequential_4_602708:
��"
sequential_4_602710:	�'
sequential_4_602712:
��"
sequential_4_602714:	�&
sequential_4_602716:	�!
sequential_4_602718:!
readvariableop_resource: #
readvariableop_1_resource: #
readvariableop_2_resource: %
sequential_5_602743:!
sequential_5_602745:
identity

identity_1

identity_2

identity_3��ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�,dense_16/bias/Regularizer/Abs/ReadVariableOp�/dense_16/bias/Regularizer/Square/ReadVariableOp�.dense_16/kernel/Regularizer/Abs/ReadVariableOp�1dense_16/kernel/Regularizer/Square/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOp�,dense_18/bias/Regularizer/Abs/ReadVariableOp�/dense_18/bias/Regularizer/Square/ReadVariableOp�.dense_18/kernel/Regularizer/Abs/ReadVariableOp�1dense_18/kernel/Regularizer/Square/ReadVariableOp�,dense_19/bias/Regularizer/Abs/ReadVariableOp�/dense_19/bias/Regularizer/Square/ReadVariableOp�.dense_19/kernel/Regularizer/Abs/ReadVariableOp�1dense_19/kernel/Regularizer/Square/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOp�,dense_21/bias/Regularizer/Abs/ReadVariableOp�/dense_21/bias/Regularizer/Square/ReadVariableOp�.dense_21/kernel/Regularizer/Abs/ReadVariableOp�1dense_21/kernel/Regularizer/Square/ReadVariableOp�,dense_22/bias/Regularizer/Abs/ReadVariableOp�/dense_22/bias/Regularizer/Square/ReadVariableOp�.dense_22/kernel/Regularizer/Abs/ReadVariableOp�1dense_22/kernel/Regularizer/Square/ReadVariableOp�,dense_23/bias/Regularizer/Abs/ReadVariableOp�/dense_23/bias/Regularizer/Square/ReadVariableOp�.dense_23/kernel/Regularizer/Abs/ReadVariableOp�1dense_23/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�&sequential_4/StatefulPartitionedCall_1�$sequential_5/StatefulPartitionedCall�&sequential_5/StatefulPartitionedCall_1�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallxsequential_4_602692sequential_4_602694sequential_4_602696sequential_4_602698sequential_4_602700sequential_4_602702sequential_4_602704sequential_4_602706sequential_4_602708sequential_4_602710sequential_4_602712sequential_4_602714sequential_4_602716sequential_4_602718*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493d
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        f
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       f
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_sliceStridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_mask^
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
: *
dtype0h
mulMulReadVariableOp:value:0strided_slice:output:0*
T0*#
_output_shapes
:���������f
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_1StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
: *
dtype0n
mul_1MulReadVariableOp_1:value:0strided_slice_1:output:0*
T0*#
_output_shapes
:���������f
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_2StridedSlice-sequential_4/StatefulPartitionedCall:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskb
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
: *
dtype0n
mul_2MulReadVariableOp_2:value:0strided_slice_2:output:0*
T0*#
_output_shapes
:���������|
stackPackmul:z:0	mul_1:z:0	mul_2:z:0*
N*
T0*'
_output_shapes
:���������*
axis����������
$sequential_5/StatefulPartitionedCallStatefulPartitionedCallstack:output:0sequential_5_602743sequential_5_602745*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524�
&sequential_5/StatefulPartitionedCall_1StatefulPartitionedCall-sequential_4/StatefulPartitionedCall:output:0sequential_5_602743sequential_5_602745*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_5_layer_call_and_return_conditional_losses_602524p
subSubx/sequential_5/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������K
SquareSquaresub:z:0*
T0*'
_output_shapes
:���������V
ConstConst*
_output_shapes
:*
dtype0*
valueB"       I
MeanMean
Square:y:0Const:output:0*
T0*
_output_shapes
: f
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_3StridedSlicexstrided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskL
mul_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?f
mul_3Mulmul_3/x:output:0strided_slice_3:output:0*
T0*#
_output_shapes
:���������f
strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_4StridedSlicexstrided_slice_4/stack:output:0 strided_slice_4/stack_1:output:0 strided_slice_4/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskf
strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        h
strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       h
strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      �
strided_slice_5StridedSlicexstrided_slice_5/stack:output:0 strided_slice_5/stack_1:output:0 strided_slice_5/stack_2:output:0*
Index0*
T0*#
_output_shapes
:���������*

begin_mask*
end_mask*
shrink_axis_maskZ
Square_1Squarestrided_slice_5:output:0*
T0*#
_output_shapes
:���������b
sub_1Substrided_slice_4:output:0Square_1:y:0*
T0*#
_output_shapes
:���������L
mul_4/xConst*
_output_shapes
: *
dtype0*
valueB
 *   ?W
mul_4Mulmul_4/x:output:0	sub_1:z:0*
T0*#
_output_shapes
:���������u
stack_1Pack	mul_3:z:0	mul_4:z:0*
N*
T0*'
_output_shapes
:���������*
axis���������
sub_2Sub-sequential_5/StatefulPartitionedCall:output:0stack_1:output:0*
T0*'
_output_shapes
:���������O
Square_2Square	sub_2:z:0*
T0*'
_output_shapes
:���������X
Const_1Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_1MeanSquare_2:y:0Const_1:output:0*
T0*
_output_shapes
: �
&sequential_4/StatefulPartitionedCall_1StatefulPartitionedCallstack_1:output:0sequential_4_602692sequential_4_602694sequential_4_602696sequential_4_602698sequential_4_602700sequential_4_602702sequential_4_602704sequential_4_602706sequential_4_602708sequential_4_602710sequential_4_602712sequential_4_602714sequential_4_602716sequential_4_602718*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601493
sub_3Substack:output:0/sequential_4/StatefulPartitionedCall_1:output:0*
T0*'
_output_shapes
:���������O
Square_3Square	sub_3:z:0*
T0*'
_output_shapes
:���������X
Const_2Const*
_output_shapes
:*
dtype0*
valueB"       O
Mean_2MeanSquare_3:y:0Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_16/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602692*
_output_shapes
:	�*
dtype0�
dense_16/kernel/Regularizer/AbsAbs6dense_16/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_16/kernel/Regularizer/SumSum#dense_16/kernel/Regularizer/Abs:y:0,dense_16/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_16/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/kernel/Regularizer/mulMul*dense_16/kernel/Regularizer/mul/x:output:0(dense_16/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/kernel/Regularizer/addAddV2*dense_16/kernel/Regularizer/Const:output:0#dense_16/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_16/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602692*
_output_shapes
:	�*
dtype0�
"dense_16/kernel/Regularizer/SquareSquare9dense_16/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_16/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_16/kernel/Regularizer/Sum_1Sum&dense_16/kernel/Regularizer/Square:y:0,dense_16/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_16/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_16/kernel/Regularizer/mul_1Mul,dense_16/kernel/Regularizer/mul_1/x:output:0*dense_16/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_16/kernel/Regularizer/add_1AddV2#dense_16/kernel/Regularizer/add:z:0%dense_16/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_16/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602694*
_output_shapes	
:�*
dtype0�
dense_16/bias/Regularizer/AbsAbs4dense_16/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/SumSum!dense_16/bias/Regularizer/Abs:y:0*dense_16/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_16/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mulMul(dense_16/bias/Regularizer/mul/x:output:0&dense_16/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/addAddV2(dense_16/bias/Regularizer/Const:output:0!dense_16/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_16/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602694*
_output_shapes	
:�*
dtype0�
 dense_16/bias/Regularizer/SquareSquare7dense_16/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_16/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_16/bias/Regularizer/Sum_1Sum$dense_16/bias/Regularizer/Square:y:0*dense_16/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_16/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_16/bias/Regularizer/mul_1Mul*dense_16/bias/Regularizer/mul_1/x:output:0(dense_16/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_16/bias/Regularizer/add_1AddV2!dense_16/bias/Regularizer/add:z:0#dense_16/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602696* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602696* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602698*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602698*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_18/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602700* 
_output_shapes
:
��*
dtype0�
dense_18/kernel/Regularizer/AbsAbs6dense_18/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_18/kernel/Regularizer/SumSum#dense_18/kernel/Regularizer/Abs:y:0,dense_18/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_18/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/kernel/Regularizer/mulMul*dense_18/kernel/Regularizer/mul/x:output:0(dense_18/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/kernel/Regularizer/addAddV2*dense_18/kernel/Regularizer/Const:output:0#dense_18/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_18/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602700* 
_output_shapes
:
��*
dtype0�
"dense_18/kernel/Regularizer/SquareSquare9dense_18/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_18/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_18/kernel/Regularizer/Sum_1Sum&dense_18/kernel/Regularizer/Square:y:0,dense_18/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_18/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_18/kernel/Regularizer/mul_1Mul,dense_18/kernel/Regularizer/mul_1/x:output:0*dense_18/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_18/kernel/Regularizer/add_1AddV2#dense_18/kernel/Regularizer/add:z:0%dense_18/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_18/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602702*
_output_shapes	
:�*
dtype0�
dense_18/bias/Regularizer/AbsAbs4dense_18/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/SumSum!dense_18/bias/Regularizer/Abs:y:0*dense_18/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_18/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mulMul(dense_18/bias/Regularizer/mul/x:output:0&dense_18/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/addAddV2(dense_18/bias/Regularizer/Const:output:0!dense_18/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_18/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602702*
_output_shapes	
:�*
dtype0�
 dense_18/bias/Regularizer/SquareSquare7dense_18/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_18/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_18/bias/Regularizer/Sum_1Sum$dense_18/bias/Regularizer/Square:y:0*dense_18/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_18/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_18/bias/Regularizer/mul_1Mul*dense_18/bias/Regularizer/mul_1/x:output:0(dense_18/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_18/bias/Regularizer/add_1AddV2!dense_18/bias/Regularizer/add:z:0#dense_18/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_19/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602704* 
_output_shapes
:
��*
dtype0�
dense_19/kernel/Regularizer/AbsAbs6dense_19/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_19/kernel/Regularizer/SumSum#dense_19/kernel/Regularizer/Abs:y:0,dense_19/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_19/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/kernel/Regularizer/mulMul*dense_19/kernel/Regularizer/mul/x:output:0(dense_19/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/kernel/Regularizer/addAddV2*dense_19/kernel/Regularizer/Const:output:0#dense_19/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_19/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602704* 
_output_shapes
:
��*
dtype0�
"dense_19/kernel/Regularizer/SquareSquare9dense_19/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_19/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_19/kernel/Regularizer/Sum_1Sum&dense_19/kernel/Regularizer/Square:y:0,dense_19/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_19/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_19/kernel/Regularizer/mul_1Mul,dense_19/kernel/Regularizer/mul_1/x:output:0*dense_19/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_19/kernel/Regularizer/add_1AddV2#dense_19/kernel/Regularizer/add:z:0%dense_19/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_19/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602706*
_output_shapes	
:�*
dtype0�
dense_19/bias/Regularizer/AbsAbs4dense_19/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/SumSum!dense_19/bias/Regularizer/Abs:y:0*dense_19/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_19/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mulMul(dense_19/bias/Regularizer/mul/x:output:0&dense_19/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/addAddV2(dense_19/bias/Regularizer/Const:output:0!dense_19/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_19/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602706*
_output_shapes	
:�*
dtype0�
 dense_19/bias/Regularizer/SquareSquare7dense_19/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_19/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_19/bias/Regularizer/Sum_1Sum$dense_19/bias/Regularizer/Square:y:0*dense_19/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_19/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_19/bias/Regularizer/mul_1Mul*dense_19/bias/Regularizer/mul_1/x:output:0(dense_19/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_19/bias/Regularizer/add_1AddV2!dense_19/bias/Regularizer/add:z:0#dense_19/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602708* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602708* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602710*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602710*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_21/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602712* 
_output_shapes
:
��*
dtype0�
dense_21/kernel/Regularizer/AbsAbs6dense_21/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_21/kernel/Regularizer/SumSum#dense_21/kernel/Regularizer/Abs:y:0,dense_21/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_21/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/kernel/Regularizer/mulMul*dense_21/kernel/Regularizer/mul/x:output:0(dense_21/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/kernel/Regularizer/addAddV2*dense_21/kernel/Regularizer/Const:output:0#dense_21/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_21/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602712* 
_output_shapes
:
��*
dtype0�
"dense_21/kernel/Regularizer/SquareSquare9dense_21/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_21/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_21/kernel/Regularizer/Sum_1Sum&dense_21/kernel/Regularizer/Square:y:0,dense_21/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_21/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_21/kernel/Regularizer/mul_1Mul,dense_21/kernel/Regularizer/mul_1/x:output:0*dense_21/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_21/kernel/Regularizer/add_1AddV2#dense_21/kernel/Regularizer/add:z:0%dense_21/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    }
,dense_21/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602714*
_output_shapes	
:�*
dtype0�
dense_21/bias/Regularizer/AbsAbs4dense_21/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/SumSum!dense_21/bias/Regularizer/Abs:y:0*dense_21/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_21/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mulMul(dense_21/bias/Regularizer/mul/x:output:0&dense_21/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/addAddV2(dense_21/bias/Regularizer/Const:output:0!dense_21/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_21/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602714*
_output_shapes	
:�*
dtype0�
 dense_21/bias/Regularizer/SquareSquare7dense_21/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_21/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_21/bias/Regularizer/Sum_1Sum$dense_21/bias/Regularizer/Square:y:0*dense_21/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_21/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_21/bias/Regularizer/mul_1Mul*dense_21/bias/Regularizer/mul_1/x:output:0(dense_21/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_21/bias/Regularizer/add_1AddV2!dense_21/bias/Regularizer/add:z:0#dense_21/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_22/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602716*
_output_shapes
:	�*
dtype0�
dense_22/kernel/Regularizer/AbsAbs6dense_22/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_22/kernel/Regularizer/SumSum#dense_22/kernel/Regularizer/Abs:y:0,dense_22/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_22/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/kernel/Regularizer/mulMul*dense_22/kernel/Regularizer/mul/x:output:0(dense_22/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/kernel/Regularizer/addAddV2*dense_22/kernel/Regularizer/Const:output:0#dense_22/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_22/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602716*
_output_shapes
:	�*
dtype0�
"dense_22/kernel/Regularizer/SquareSquare9dense_22/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�t
#dense_22/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_22/kernel/Regularizer/Sum_1Sum&dense_22/kernel/Regularizer/Square:y:0,dense_22/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_22/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_22/kernel/Regularizer/mul_1Mul,dense_22/kernel/Regularizer/mul_1/x:output:0*dense_22/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_22/kernel/Regularizer/add_1AddV2#dense_22/kernel/Regularizer/add:z:0%dense_22/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_22/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_4_602718*
_output_shapes
:*
dtype0
dense_22/bias/Regularizer/AbsAbs4dense_22/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/SumSum!dense_22/bias/Regularizer/Abs:y:0*dense_22/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_22/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mulMul(dense_22/bias/Regularizer/mul/x:output:0&dense_22/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/addAddV2(dense_22/bias/Regularizer/Const:output:0!dense_22/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_22/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_602718*
_output_shapes
:*
dtype0�
 dense_22/bias/Regularizer/SquareSquare7dense_22/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_22/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_22/bias/Regularizer/Sum_1Sum$dense_22/bias/Regularizer/Square:y:0*dense_22/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_22/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_22/bias/Regularizer/mul_1Mul*dense_22/bias/Regularizer/mul_1/x:output:0(dense_22/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_22/bias/Regularizer/add_1AddV2!dense_22/bias/Regularizer/add:z:0#dense_22/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_23/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_602743*
_output_shapes

:*
dtype0�
dense_23/kernel/Regularizer/AbsAbs6dense_23/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_23/kernel/Regularizer/SumSum#dense_23/kernel/Regularizer/Abs:y:0,dense_23/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_23/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/kernel/Regularizer/mulMul*dense_23/kernel/Regularizer/mul/x:output:0(dense_23/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/kernel/Regularizer/addAddV2*dense_23/kernel/Regularizer/Const:output:0#dense_23/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_23/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_602743*
_output_shapes

:*
dtype0�
"dense_23/kernel/Regularizer/SquareSquare9dense_23/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:t
#dense_23/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_23/kernel/Regularizer/Sum_1Sum&dense_23/kernel/Regularizer/Square:y:0,dense_23/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_23/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_23/kernel/Regularizer/mul_1Mul,dense_23/kernel/Regularizer/mul_1/x:output:0*dense_23/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_23/kernel/Regularizer/add_1AddV2#dense_23/kernel/Regularizer/add:z:0%dense_23/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    |
,dense_23/bias/Regularizer/Abs/ReadVariableOpReadVariableOpsequential_5_602745*
_output_shapes
:*
dtype0
dense_23/bias/Regularizer/AbsAbs4dense_23/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/SumSum!dense_23/bias/Regularizer/Abs:y:0*dense_23/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_23/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mulMul(dense_23/bias/Regularizer/mul/x:output:0&dense_23/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/addAddV2(dense_23/bias/Regularizer/Const:output:0!dense_23/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: 
/dense_23/bias/Regularizer/Square/ReadVariableOpReadVariableOpsequential_5_602745*
_output_shapes
:*
dtype0�
 dense_23/bias/Regularizer/SquareSquare7dense_23/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:k
!dense_23/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_23/bias/Regularizer/Sum_1Sum$dense_23/bias/Regularizer/Square:y:0*dense_23/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_23/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_23/bias/Regularizer/mul_1Mul*dense_23/bias/Regularizer/mul_1/x:output:0(dense_23/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_23/bias/Regularizer/add_1AddV2!dense_23/bias/Regularizer/add:z:0#dense_23/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: |
IdentityIdentity-sequential_5/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������M

Identity_1IdentityMean:output:0^NoOp*
T0*
_output_shapes
: O

Identity_2IdentityMean_1:output:0^NoOp*
T0*
_output_shapes
: O

Identity_3IdentityMean_2:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^ReadVariableOp^ReadVariableOp_1^ReadVariableOp_2-^dense_16/bias/Regularizer/Abs/ReadVariableOp0^dense_16/bias/Regularizer/Square/ReadVariableOp/^dense_16/kernel/Regularizer/Abs/ReadVariableOp2^dense_16/kernel/Regularizer/Square/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp-^dense_18/bias/Regularizer/Abs/ReadVariableOp0^dense_18/bias/Regularizer/Square/ReadVariableOp/^dense_18/kernel/Regularizer/Abs/ReadVariableOp2^dense_18/kernel/Regularizer/Square/ReadVariableOp-^dense_19/bias/Regularizer/Abs/ReadVariableOp0^dense_19/bias/Regularizer/Square/ReadVariableOp/^dense_19/kernel/Regularizer/Abs/ReadVariableOp2^dense_19/kernel/Regularizer/Square/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp-^dense_21/bias/Regularizer/Abs/ReadVariableOp0^dense_21/bias/Regularizer/Square/ReadVariableOp/^dense_21/kernel/Regularizer/Abs/ReadVariableOp2^dense_21/kernel/Regularizer/Square/ReadVariableOp-^dense_22/bias/Regularizer/Abs/ReadVariableOp0^dense_22/bias/Regularizer/Square/ReadVariableOp/^dense_22/kernel/Regularizer/Abs/ReadVariableOp2^dense_22/kernel/Regularizer/Square/ReadVariableOp-^dense_23/bias/Regularizer/Abs/ReadVariableOp0^dense_23/bias/Regularizer/Square/ReadVariableOp/^dense_23/kernel/Regularizer/Abs/ReadVariableOp2^dense_23/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall'^sequential_4/StatefulPartitionedCall_1%^sequential_5/StatefulPartitionedCall'^sequential_5/StatefulPartitionedCall_1*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*L
_input_shapes;
9:���������: : : : : : : : : : : : : : : : : : : 2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22\
,dense_16/bias/Regularizer/Abs/ReadVariableOp,dense_16/bias/Regularizer/Abs/ReadVariableOp2b
/dense_16/bias/Regularizer/Square/ReadVariableOp/dense_16/bias/Regularizer/Square/ReadVariableOp2`
.dense_16/kernel/Regularizer/Abs/ReadVariableOp.dense_16/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_16/kernel/Regularizer/Square/ReadVariableOp1dense_16/kernel/Regularizer/Square/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp2\
,dense_18/bias/Regularizer/Abs/ReadVariableOp,dense_18/bias/Regularizer/Abs/ReadVariableOp2b
/dense_18/bias/Regularizer/Square/ReadVariableOp/dense_18/bias/Regularizer/Square/ReadVariableOp2`
.dense_18/kernel/Regularizer/Abs/ReadVariableOp.dense_18/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_18/kernel/Regularizer/Square/ReadVariableOp1dense_18/kernel/Regularizer/Square/ReadVariableOp2\
,dense_19/bias/Regularizer/Abs/ReadVariableOp,dense_19/bias/Regularizer/Abs/ReadVariableOp2b
/dense_19/bias/Regularizer/Square/ReadVariableOp/dense_19/bias/Regularizer/Square/ReadVariableOp2`
.dense_19/kernel/Regularizer/Abs/ReadVariableOp.dense_19/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_19/kernel/Regularizer/Square/ReadVariableOp1dense_19/kernel/Regularizer/Square/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp2\
,dense_21/bias/Regularizer/Abs/ReadVariableOp,dense_21/bias/Regularizer/Abs/ReadVariableOp2b
/dense_21/bias/Regularizer/Square/ReadVariableOp/dense_21/bias/Regularizer/Square/ReadVariableOp2`
.dense_21/kernel/Regularizer/Abs/ReadVariableOp.dense_21/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_21/kernel/Regularizer/Square/ReadVariableOp1dense_21/kernel/Regularizer/Square/ReadVariableOp2\
,dense_22/bias/Regularizer/Abs/ReadVariableOp,dense_22/bias/Regularizer/Abs/ReadVariableOp2b
/dense_22/bias/Regularizer/Square/ReadVariableOp/dense_22/bias/Regularizer/Square/ReadVariableOp2`
.dense_22/kernel/Regularizer/Abs/ReadVariableOp.dense_22/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_22/kernel/Regularizer/Square/ReadVariableOp1dense_22/kernel/Regularizer/Square/ReadVariableOp2\
,dense_23/bias/Regularizer/Abs/ReadVariableOp,dense_23/bias/Regularizer/Abs/ReadVariableOp2b
/dense_23/bias/Regularizer/Square/ReadVariableOp/dense_23/bias/Regularizer/Square/ReadVariableOp2`
.dense_23/kernel/Regularizer/Abs/ReadVariableOp.dense_23/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_23/kernel/Regularizer/Square/ReadVariableOp1dense_23/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall2P
&sequential_4/StatefulPartitionedCall_1&sequential_4/StatefulPartitionedCall_12L
$sequential_5/StatefulPartitionedCall$sequential_5/StatefulPartitionedCall2P
&sequential_5/StatefulPartitionedCall_1&sequential_5/StatefulPartitionedCall_1:J F
'
_output_shapes
:���������

_user_specified_namex
�0
�
D__inference_dense_20_layer_call_and_return_conditional_losses_601182

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOp�.dense_20/kernel/Regularizer/Abs/ReadVariableOp�1dense_20/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_20/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_20/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_20/kernel/Regularizer/AbsAbs6dense_20/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_20/kernel/Regularizer/SumSum#dense_20/kernel/Regularizer/Abs:y:0,dense_20/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_20/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/kernel/Regularizer/mulMul*dense_20/kernel/Regularizer/mul/x:output:0(dense_20/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/kernel/Regularizer/addAddV2*dense_20/kernel/Regularizer/Const:output:0#dense_20/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_20/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_20/kernel/Regularizer/SquareSquare9dense_20/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_20/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_20/kernel/Regularizer/Sum_1Sum&dense_20/kernel/Regularizer/Square:y:0,dense_20/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_20/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_20/kernel/Regularizer/mul_1Mul,dense_20/kernel/Regularizer/mul_1/x:output:0*dense_20/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_20/kernel/Regularizer/add_1AddV2#dense_20/kernel/Regularizer/add:z:0%dense_20/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp/^dense_20/kernel/Regularizer/Abs/ReadVariableOp2^dense_20/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp2`
.dense_20/kernel/Regularizer/Abs/ReadVariableOp.dense_20/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_20/kernel/Regularizer/Square/ReadVariableOp1dense_20/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
)__inference_dense_22_layer_call_fn_606899

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dense_22_layer_call_and_return_conditional_losses_601276o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
__inference_loss_fn_9_607140D
5dense_20_bias_regularizer_abs_readvariableop_resource:	�
identity��,dense_20/bias/Regularizer/Abs/ReadVariableOp�/dense_20/bias/Regularizer/Square/ReadVariableOpd
dense_20/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_20/bias/Regularizer/Abs/ReadVariableOpReadVariableOp5dense_20_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_20/bias/Regularizer/AbsAbs4dense_20/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/SumSum!dense_20/bias/Regularizer/Abs:y:0*dense_20/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_20/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mulMul(dense_20/bias/Regularizer/mul/x:output:0&dense_20/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/addAddV2(dense_20/bias/Regularizer/Const:output:0!dense_20/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_20/bias/Regularizer/Square/ReadVariableOpReadVariableOp5dense_20_bias_regularizer_abs_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_20/bias/Regularizer/SquareSquare7dense_20/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_20/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_20/bias/Regularizer/Sum_1Sum$dense_20/bias/Regularizer/Square:y:0*dense_20/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_20/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_20/bias/Regularizer/mul_1Mul*dense_20/bias/Regularizer/mul_1/x:output:0(dense_20/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_20/bias/Regularizer/add_1AddV2!dense_20/bias/Regularizer/add:z:0#dense_20/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: a
IdentityIdentity#dense_20/bias/Regularizer/add_1:z:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp-^dense_20/bias/Regularizer/Abs/ReadVariableOp0^dense_20/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2\
,dense_20/bias/Regularizer/Abs/ReadVariableOp,dense_20/bias/Regularizer/Abs/ReadVariableOp2b
/dense_20/bias/Regularizer/Square/ReadVariableOp/dense_20/bias/Regularizer/Square/ReadVariableOp
�
�
-__inference_sequential_4_layer_call_fn_605726

inputs
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:
��
	unknown_4:	�
	unknown_5:
��
	unknown_6:	�
	unknown_7:
��
	unknown_8:	�
	unknown_9:
��

unknown_10:	�

unknown_11:	�

unknown_12:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*0
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *Q
fLRJ
H__inference_sequential_4_layer_call_and_return_conditional_losses_601878o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*B
_input_shapes1
/:���������: : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�0
�
D__inference_dense_17_layer_call_and_return_conditional_losses_606540

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�,dense_17/bias/Regularizer/Abs/ReadVariableOp�/dense_17/bias/Regularizer/Square/ReadVariableOp�.dense_17/kernel/Regularizer/Abs/ReadVariableOp�1dense_17/kernel/Regularizer/Square/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
SeluSeluBiasAdd:output:0*
T0*(
_output_shapes
:����������f
!dense_17/kernel/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
.dense_17/kernel/Regularizer/Abs/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
dense_17/kernel/Regularizer/AbsAbs6dense_17/kernel/Regularizer/Abs/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB"       �
dense_17/kernel/Regularizer/SumSum#dense_17/kernel/Regularizer/Abs:y:0,dense_17/kernel/Regularizer/Const_1:output:0*
T0*
_output_shapes
: f
!dense_17/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/kernel/Regularizer/mulMul*dense_17/kernel/Regularizer/mul/x:output:0(dense_17/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/kernel/Regularizer/addAddV2*dense_17/kernel/Regularizer/Const:output:0#dense_17/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: �
1dense_17/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
"dense_17/kernel/Regularizer/SquareSquare9dense_17/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��t
#dense_17/kernel/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB"       �
!dense_17/kernel/Regularizer/Sum_1Sum&dense_17/kernel/Regularizer/Square:y:0,dense_17/kernel/Regularizer/Const_2:output:0*
T0*
_output_shapes
: h
#dense_17/kernel/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
!dense_17/kernel/Regularizer/mul_1Mul,dense_17/kernel/Regularizer/mul_1/x:output:0*dense_17/kernel/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
!dense_17/kernel/Regularizer/add_1AddV2#dense_17/kernel/Regularizer/add:z:0%dense_17/kernel/Regularizer/mul_1:z:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    �
,dense_17/bias/Regularizer/Abs/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
dense_17/bias/Regularizer/AbsAbs4dense_17/bias/Regularizer/Abs/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_1Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/SumSum!dense_17/bias/Regularizer/Abs:y:0*dense_17/bias/Regularizer/Const_1:output:0*
T0*
_output_shapes
: d
dense_17/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mulMul(dense_17/bias/Regularizer/mul/x:output:0&dense_17/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/addAddV2(dense_17/bias/Regularizer/Const:output:0!dense_17/bias/Regularizer/mul:z:0*
T0*
_output_shapes
: �
/dense_17/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 dense_17/bias/Regularizer/SquareSquare7dense_17/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:�k
!dense_17/bias/Regularizer/Const_2Const*
_output_shapes
:*
dtype0*
valueB: �
dense_17/bias/Regularizer/Sum_1Sum$dense_17/bias/Regularizer/Square:y:0*dense_17/bias/Regularizer/Const_2:output:0*
T0*
_output_shapes
: f
!dense_17/bias/Regularizer/mul_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *}�&�
dense_17/bias/Regularizer/mul_1Mul*dense_17/bias/Regularizer/mul_1/x:output:0(dense_17/bias/Regularizer/Sum_1:output:0*
T0*
_output_shapes
: �
dense_17/bias/Regularizer/add_1AddV2!dense_17/bias/Regularizer/add:z:0#dense_17/bias/Regularizer/mul_1:z:0*
T0*
_output_shapes
: b
IdentityIdentitySelu:activations:0^NoOp*
T0*(
_output_shapes
:�����������
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp-^dense_17/bias/Regularizer/Abs/ReadVariableOp0^dense_17/bias/Regularizer/Square/ReadVariableOp/^dense_17/kernel/Regularizer/Abs/ReadVariableOp2^dense_17/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2\
,dense_17/bias/Regularizer/Abs/ReadVariableOp,dense_17/bias/Regularizer/Abs/ReadVariableOp2b
/dense_17/bias/Regularizer/Square/ReadVariableOp/dense_17/bias/Regularizer/Square/ReadVariableOp2`
.dense_17/kernel/Regularizer/Abs/ReadVariableOp.dense_17/kernel/Regularizer/Abs/ReadVariableOp2f
1dense_17/kernel/Regularizer/Square/ReadVariableOp1dense_17/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
;
input_10
serving_default_input_1:0���������<
output_10
StatefulPartitionedCall:0���������tensorflow/serving/predict:��
�
a1
a2
a3
encoder
decoder
	optimizer
	variables
trainable_variables
	regularization_losses

	keras_api

signatures
�__call__
+�&call_and_return_all_conditional_losses
�_default_save_signature"
_tf_keras_model
: 2Variable
: 2Variable
: 2Variable
�
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
layer_with_weights-4
layer-4
layer_with_weights-5
layer-5
layer_with_weights-6
layer-6
	variables
trainable_variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_sequential
�
layer_with_weights-0
layer-0
	variables
trainable_variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_sequential
�
iter

beta_1

beta_2
	decay
 learning_rate!m�"m�#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�-m�.m�/m�0m�!v�"v�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�-v�.v�/v�0v�"
	optimizer
�
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015
16
17
18"
trackable_list_wrapper
�
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13
/14
015"
trackable_list_wrapper
 "
trackable_list_wrapper
�
1non_trainable_variables

2layers
3metrics
4layer_regularization_losses
5layer_metrics
	variables
trainable_variables
	regularization_losses
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
�

!kernel
"bias
6	variables
7trainable_variables
8regularization_losses
9	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

#kernel
$bias
:	variables
;trainable_variables
<regularization_losses
=	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

%kernel
&bias
>	variables
?trainable_variables
@regularization_losses
A	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

'kernel
(bias
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

)kernel
*bias
F	variables
Gtrainable_variables
Hregularization_losses
I	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

+kernel
,bias
J	variables
Ktrainable_variables
Lregularization_losses
M	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�

-kernel
.bias
N	variables
Otrainable_variables
Pregularization_losses
Q	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
�
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13"
trackable_list_wrapper
�
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
-12
.13"
trackable_list_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13"
trackable_list_wrapper
�
Rnon_trainable_variables

Slayers
Tmetrics
Ulayer_regularization_losses
Vlayer_metrics
	variables
trainable_variables
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�

/kernel
0bias
W	variables
Xtrainable_variables
Yregularization_losses
Z	keras_api
�__call__
+�&call_and_return_all_conditional_losses"
_tf_keras_layer
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
[non_trainable_variables

\layers
]metrics
^layer_regularization_losses
_layer_metrics
	variables
trainable_variables
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
": 	�2dense_16/kernel
:�2dense_16/bias
#:!
��2dense_17/kernel
:�2dense_17/bias
#:!
��2dense_18/kernel
:�2dense_18/bias
#:!
��2dense_19/kernel
:�2dense_19/bias
#:!
��2dense_20/kernel
:�2dense_20/bias
#:!
��2dense_21/kernel
:�2dense_21/bias
": 	�2dense_22/kernel
:2dense_22/bias
!:2dense_23/kernel
:2dense_23/bias
5
0
1
2"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
'
`0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
!0
"1"
trackable_list_wrapper
.
!0
"1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
anon_trainable_variables

blayers
cmetrics
dlayer_regularization_losses
elayer_metrics
6	variables
7trainable_variables
8regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
#0
$1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
fnon_trainable_variables

glayers
hmetrics
ilayer_regularization_losses
jlayer_metrics
:	variables
;trainable_variables
<regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
%0
&1"
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
knon_trainable_variables

llayers
mmetrics
nlayer_regularization_losses
olayer_metrics
>	variables
?trainable_variables
@regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
'0
(1"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
pnon_trainable_variables

qlayers
rmetrics
slayer_regularization_losses
tlayer_metrics
B	variables
Ctrainable_variables
Dregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
unon_trainable_variables

vlayers
wmetrics
xlayer_regularization_losses
ylayer_metrics
F	variables
Gtrainable_variables
Hregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
znon_trainable_variables

{layers
|metrics
}layer_regularization_losses
~layer_metrics
J	variables
Ktrainable_variables
Lregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
.
-0
.1"
trackable_list_wrapper
.
-0
.1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
N	variables
Otrainable_variables
Pregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
Q
0
1
2
3
4
5
6"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
W	variables
Xtrainable_variables
Yregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
R

�total

�count
�	variables
�	keras_api"
_tf_keras_metric
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_dict_wrapper
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
':%	�2Adam/dense_16/kernel/m
!:�2Adam/dense_16/bias/m
(:&
��2Adam/dense_17/kernel/m
!:�2Adam/dense_17/bias/m
(:&
��2Adam/dense_18/kernel/m
!:�2Adam/dense_18/bias/m
(:&
��2Adam/dense_19/kernel/m
!:�2Adam/dense_19/bias/m
(:&
��2Adam/dense_20/kernel/m
!:�2Adam/dense_20/bias/m
(:&
��2Adam/dense_21/kernel/m
!:�2Adam/dense_21/bias/m
':%	�2Adam/dense_22/kernel/m
 :2Adam/dense_22/bias/m
&:$2Adam/dense_23/kernel/m
 :2Adam/dense_23/bias/m
':%	�2Adam/dense_16/kernel/v
!:�2Adam/dense_16/bias/v
(:&
��2Adam/dense_17/kernel/v
!:�2Adam/dense_17/bias/v
(:&
��2Adam/dense_18/kernel/v
!:�2Adam/dense_18/bias/v
(:&
��2Adam/dense_19/kernel/v
!:�2Adam/dense_19/bias/v
(:&
��2Adam/dense_20/kernel/v
!:�2Adam/dense_20/bias/v
(:&
��2Adam/dense_21/kernel/v
!:�2Adam/dense_21/bias/v
':%	�2Adam/dense_22/kernel/v
 :2Adam/dense_22/bias/v
&:$2Adam/dense_23/kernel/v
 :2Adam/dense_23/bias/v
�2�
,__inference_conjugacy_2_layer_call_fn_603085
,__inference_conjugacy_2_layer_call_fn_604616
,__inference_conjugacy_2_layer_call_fn_604662
,__inference_conjugacy_2_layer_call_fn_603575�
���
FullArgSpec$
args�
jself
jx

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605056
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605450
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603927
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_604279�
���
FullArgSpec$
args�
jself
jx

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
!__inference__wrapped_model_600946input_1"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
-__inference_sequential_4_layer_call_fn_601524
-__inference_sequential_4_layer_call_fn_605693
-__inference_sequential_4_layer_call_fn_605726
-__inference_sequential_4_layer_call_fn_601942�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
H__inference_sequential_4_layer_call_and_return_conditional_losses_605989
H__inference_sequential_4_layer_call_and_return_conditional_losses_606252
H__inference_sequential_4_layer_call_and_return_conditional_losses_602191
H__inference_sequential_4_layer_call_and_return_conditional_losses_602440�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
-__inference_sequential_5_layer_call_fn_602531
-__inference_sequential_5_layer_call_fn_606291
-__inference_sequential_5_layer_call_fn_606300
-__inference_sequential_5_layer_call_fn_602607�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
H__inference_sequential_5_layer_call_and_return_conditional_losses_606340
H__inference_sequential_5_layer_call_and_return_conditional_losses_606380
H__inference_sequential_5_layer_call_and_return_conditional_losses_602646
H__inference_sequential_5_layer_call_and_return_conditional_losses_602685�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�B�
$__inference_signature_wrapper_604570input_1"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_16_layer_call_fn_606419�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_16_layer_call_and_return_conditional_losses_606460�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_17_layer_call_fn_606499�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_17_layer_call_and_return_conditional_losses_606540�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_18_layer_call_fn_606579�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_18_layer_call_and_return_conditional_losses_606620�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_19_layer_call_fn_606659�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_19_layer_call_and_return_conditional_losses_606700�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_20_layer_call_fn_606739�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_20_layer_call_and_return_conditional_losses_606780�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_21_layer_call_fn_606819�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_21_layer_call_and_return_conditional_losses_606860�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dense_22_layer_call_fn_606899�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_22_layer_call_and_return_conditional_losses_606940�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
__inference_loss_fn_0_606960�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_1_606980�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_2_607000�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_3_607020�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_4_607040�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_5_607060�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_6_607080�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_7_607100�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_8_607120�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_9_607140�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_10_607160�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_11_607180�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_12_607200�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_13_607220�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
)__inference_dense_23_layer_call_fn_607259�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_dense_23_layer_call_and_return_conditional_losses_607299�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
__inference_loss_fn_14_607319�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_15_607339�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� �
!__inference__wrapped_model_600946|!"#$%&'()*+,-./00�-
&�#
!�
input_1���������
� "3�0
.
output_1"�
output_1����������
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_603927�!"#$%&'()*+,-./04�1
*�'
!�
input_1���������
p 
� "O�L
�
0���������
-�*
�	
1/0 
�	
1/1 
�	
1/2 �
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_604279�!"#$%&'()*+,-./04�1
*�'
!�
input_1���������
p
� "O�L
�
0���������
-�*
�	
1/0 
�	
1/1 
�	
1/2 �
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605056�!"#$%&'()*+,-./0.�+
$�!
�
x���������
p 
� "O�L
�
0���������
-�*
�	
1/0 
�	
1/1 
�	
1/2 �
G__inference_conjugacy_2_layer_call_and_return_conditional_losses_605450�!"#$%&'()*+,-./0.�+
$�!
�
x���������
p
� "O�L
�
0���������
-�*
�	
1/0 
�	
1/1 
�	
1/2 �
,__inference_conjugacy_2_layer_call_fn_603085e!"#$%&'()*+,-./04�1
*�'
!�
input_1���������
p 
� "�����������
,__inference_conjugacy_2_layer_call_fn_603575e!"#$%&'()*+,-./04�1
*�'
!�
input_1���������
p
� "�����������
,__inference_conjugacy_2_layer_call_fn_604616_!"#$%&'()*+,-./0.�+
$�!
�
x���������
p 
� "�����������
,__inference_conjugacy_2_layer_call_fn_604662_!"#$%&'()*+,-./0.�+
$�!
�
x���������
p
� "�����������
D__inference_dense_16_layer_call_and_return_conditional_losses_606460]!"/�,
%�"
 �
inputs���������
� "&�#
�
0����������
� }
)__inference_dense_16_layer_call_fn_606419P!"/�,
%�"
 �
inputs���������
� "������������
D__inference_dense_17_layer_call_and_return_conditional_losses_606540^#$0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� ~
)__inference_dense_17_layer_call_fn_606499Q#$0�-
&�#
!�
inputs����������
� "������������
D__inference_dense_18_layer_call_and_return_conditional_losses_606620^%&0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� ~
)__inference_dense_18_layer_call_fn_606579Q%&0�-
&�#
!�
inputs����������
� "������������
D__inference_dense_19_layer_call_and_return_conditional_losses_606700^'(0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� ~
)__inference_dense_19_layer_call_fn_606659Q'(0�-
&�#
!�
inputs����������
� "������������
D__inference_dense_20_layer_call_and_return_conditional_losses_606780^)*0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� ~
)__inference_dense_20_layer_call_fn_606739Q)*0�-
&�#
!�
inputs����������
� "������������
D__inference_dense_21_layer_call_and_return_conditional_losses_606860^+,0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� ~
)__inference_dense_21_layer_call_fn_606819Q+,0�-
&�#
!�
inputs����������
� "������������
D__inference_dense_22_layer_call_and_return_conditional_losses_606940]-.0�-
&�#
!�
inputs����������
� "%�"
�
0���������
� }
)__inference_dense_22_layer_call_fn_606899P-.0�-
&�#
!�
inputs����������
� "�����������
D__inference_dense_23_layer_call_and_return_conditional_losses_607299\/0/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� |
)__inference_dense_23_layer_call_fn_607259O/0/�,
%�"
 �
inputs���������
� "����������;
__inference_loss_fn_0_606960!�

� 
� "� <
__inference_loss_fn_10_607160+�

� 
� "� <
__inference_loss_fn_11_607180,�

� 
� "� <
__inference_loss_fn_12_607200-�

� 
� "� <
__inference_loss_fn_13_607220.�

� 
� "� <
__inference_loss_fn_14_607319/�

� 
� "� <
__inference_loss_fn_15_6073390�

� 
� "� ;
__inference_loss_fn_1_606980"�

� 
� "� ;
__inference_loss_fn_2_607000#�

� 
� "� ;
__inference_loss_fn_3_607020$�

� 
� "� ;
__inference_loss_fn_4_607040%�

� 
� "� ;
__inference_loss_fn_5_607060&�

� 
� "� ;
__inference_loss_fn_6_607080'�

� 
� "� ;
__inference_loss_fn_7_607100(�

� 
� "� ;
__inference_loss_fn_8_607120)�

� 
� "� ;
__inference_loss_fn_9_607140*�

� 
� "� �
H__inference_sequential_4_layer_call_and_return_conditional_losses_602191x!"#$%&'()*+,-.?�<
5�2
(�%
dense_16_input���������
p 

 
� "%�"
�
0���������
� �
H__inference_sequential_4_layer_call_and_return_conditional_losses_602440x!"#$%&'()*+,-.?�<
5�2
(�%
dense_16_input���������
p

 
� "%�"
�
0���������
� �
H__inference_sequential_4_layer_call_and_return_conditional_losses_605989p!"#$%&'()*+,-.7�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
H__inference_sequential_4_layer_call_and_return_conditional_losses_606252p!"#$%&'()*+,-.7�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
-__inference_sequential_4_layer_call_fn_601524k!"#$%&'()*+,-.?�<
5�2
(�%
dense_16_input���������
p 

 
� "�����������
-__inference_sequential_4_layer_call_fn_601942k!"#$%&'()*+,-.?�<
5�2
(�%
dense_16_input���������
p

 
� "�����������
-__inference_sequential_4_layer_call_fn_605693c!"#$%&'()*+,-.7�4
-�*
 �
inputs���������
p 

 
� "�����������
-__inference_sequential_4_layer_call_fn_605726c!"#$%&'()*+,-.7�4
-�*
 �
inputs���������
p

 
� "�����������
H__inference_sequential_5_layer_call_and_return_conditional_losses_602646l/0?�<
5�2
(�%
dense_23_input���������
p 

 
� "%�"
�
0���������
� �
H__inference_sequential_5_layer_call_and_return_conditional_losses_602685l/0?�<
5�2
(�%
dense_23_input���������
p

 
� "%�"
�
0���������
� �
H__inference_sequential_5_layer_call_and_return_conditional_losses_606340d/07�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
H__inference_sequential_5_layer_call_and_return_conditional_losses_606380d/07�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
-__inference_sequential_5_layer_call_fn_602531_/0?�<
5�2
(�%
dense_23_input���������
p 

 
� "�����������
-__inference_sequential_5_layer_call_fn_602607_/0?�<
5�2
(�%
dense_23_input���������
p

 
� "�����������
-__inference_sequential_5_layer_call_fn_606291W/07�4
-�*
 �
inputs���������
p 

 
� "�����������
-__inference_sequential_5_layer_call_fn_606300W/07�4
-�*
 �
inputs���������
p

 
� "�����������
$__inference_signature_wrapper_604570�!"#$%&'()*+,-./0;�8
� 
1�.
,
input_1!�
input_1���������"3�0
.
output_1"�
output_1���������