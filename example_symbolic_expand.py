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

# Constructs the linear system in Nister's 5pt algorithm
#
# Works with Sage Version 5.2


from sage.all import *
from numpy.core.fromnumeric import var
from utilsym import *


def simplifyExpanded( input ):
    """Reduces number of multiplications by moving the most common elements outside.
       Example: 'a*b*c + a*d*e' will become 'a*(b*c + d*e)
    """

    if not len(input): return ''

    def removeFirst( text , var ):
        l =len(var)
        i = text.find(var)
        if i-1>0 and text[i-1] == '*':
            i -= 1;l += 1
        if i+l<len(text) and text[i+l] == '*':
            l += 1
        if l==5: l=4

        return text[0:i] + text[i+l:]

    def reconstruct( sequence ):
        output = sequence[0]
        for w in sequence[1:]:
            if w[0] == '-':
                output += ' - '+w[1:]
            else:
                output += ' + '+w
        return output

    s = input.replace('- ','-').split(' ')
    s = [w for w in s if len(w) > 1 ]

    # Find the frequency of each variable
    dict = {}
    for w in s:
        vars = w.replace('-','').split('*')
        for v in vars:
            if dict.has_key(v):
                dict[v] += 1
            else: dict[v] = 1

    bestVar = ''
    bestCount = 0
    for k,v in dict.items():
        if v > bestCount:
            bestCount = v
            bestVar = k

    include = []
    exclude = []

    for w in s:
        if bestVar in w:
            include.append(w)
        else:
            exclude.append(w)

    include = [removeFirst(w,bestVar) for w in include]

    output = bestVar + '*( '
    if bestVar[0].isdigit(): output += simplifyExpanded(reconstruct(include))
    else: output += reconstruct(include)
    output += ' )'

    if len(exclude):
        output += ' + '+simplifyExpanded( reconstruct(exclude))

    return output

def printData( var , eqs , keys ):
    """Prints a java equation"""
    f = open('%s.txt'%var,'w')
    for row,eq in enumerate(eqs):
        index = len(keys)*row
        for k in keys:
            f.write('%s.data[%d] = %s;\n'%(var,index,simplifyExpanded(extractVarEq(eq,k))))
            index += 1
    f.close()

x,y,z = var('x','y','z')

X = symMatrix( 3, 3 , 'X')
Y = symMatrix( 3, 3 , 'Y')
Z = symMatrix( 3, 3 , 'Z')
W = symMatrix( 3, 3 , 'W')

E = x*X + y*Y + z*Z + 1*W
EE = E*E.T

eq1=det(E)
eq2=EE*E-0.5*EE.trace()*E

eqs = (eq1,eq2[0,0],eq2[0,1],eq2[0,2],eq2[1,0],eq2[1,1],eq2[1,2],eq2[2,0],eq2[2,1],eq2[2,2])

keysA = ('x^3','y^3','x^2*y','x*y^2','x^2*z','x^2','y^2*z','y^2','x*y*z','x*y')
keysB = ('x*z^2','x*z','x','y*z^2','y*z','y','z^3','z^2','z','')

# print out machine code for the linear system
printData('A',eqs,keysA)
printData('B',eqs,keysB)
