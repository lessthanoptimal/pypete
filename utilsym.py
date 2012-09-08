# Copyright (c) 2012, Peter Abeles. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Utility functions for handling symbolic math and converting into Java code
#
# Works with Sage Version 5.2

from numpy.core.fromnumeric import var
from numpy.linalg.linalg import det

from sage.all import *

def symMatrix( numRows , numCols , letter ):
  """Creates a matrix where each element is a unique symbol
  """

  A = matrix(SR,numRows,numCols)
  for i in range(0,numRows):
    for j in range(0,numCols):
      A[i,j] = var('%s%d%d' % (letter,i,j) )

  return A

def expandPower( expression ):
  """Expands out variables which are multiplied by an integer value
  For example x^2 = x*x and x^4 = x*x*x*x
  """

  l = expression.split('*')

  # handle special case with a minus sign in front
  if l[0][0] == '-':
    l[0] = l[0][1:]
    expression = '-'
  else:
    expression = ''

  for s in l:
    if len(s) >= 3 and s[-2] == '^':
       var = s[:-2]
       expanded = var
       for i in range(int(s[-1])-1): expanded += '*'+var
       expression += expanded + '*'
    else:
       expression += s + '*'
  return expression[:-1]

def extractVarEq( expression , key ):
  """Expands the expression out and searches for all blocks of multiplication that are multiplied by the key.
  All other blocks are discarded and the key is removed from the selected blocks
  Example:  "x*a*b + y*a*a*b - c*x*d"  would output "a*b - c*d" if the key was 'x'
  """

  chars = set('xyz')
  # expand out and convert into a string
  expression = str(expression.expand())
  # Make sure negative symbols are not stripped and split into multiplicative blocks
  s = expression.replace('- ','-').split(' ')
  # Find blocks multiplied by the key and remove the key from the string
  if len(key) == 0:
    var = [w for w in s if len(w) != 1 and not any( c in chars for c in w )]
  else:
    var = [w[:-(1+len(key))] for w in s if (w.endswith(key) and not any(c in chars for c in w[:-len(key)])) ]
 
  # Expand out power
  var = [expandPower(w) for w in var] 
 
  # construct a string which can be compiled
  ret = var[0]
  for w in var[1:]:
    if w[0] == '-':
      ret += ' - '+w[1:]
    else:
      ret += ' + '+w

  return ret
      



