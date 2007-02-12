import string
import periodictable

def listnucleus(keystring='mass'):
    items = periodictable.isotopes.items()
    x=[[value[keystring],key] for key,value in items]
    x.sort()
    print '<table width="80%" border=1>'
    
    print '<tr>',
    print '<th align="center">Atom/isotope</th>',
    print '<th align="center">Spin</th>',
    print '<th align="center">Larmor frequency<br>(MHz)<br>at 9.396T</th>',
    print '<th align="center">Abundance<br>(%)</th>',
    print '<th align="center">Absolute<br>sensitivity</th>',
    print '<th align="center">Relative<br>sensitivity</th>',
    print '<th align="center">Quadrupole moment<br>(10<sup>-28</sup>m<sup>2</sup>)</th>',
    print '<th align="center">Ref</th>'
    print '</tr>'
    print ' '
    for key in x:
       n=periodictable.isotopes[key[1]]
       print '<tr>',
       print '<td align="center">','<span id="'+key[1]+'"></span>',n["name"],'</td>',
       print '<td align="center">',n["spin"],'</td>',
       print '<td align="center">',n["larmor"]*4,'</td>',
       print '<td align="center">',n["abundance"],'</td>',
       print '<td align="center">',n["absolute"],'</td>',
       print '<td align="center">',n["relative"],'</td>',
       print '<td align="center">',n["quadrupole"],'</td>',
       name=string.lower(n["name"].split()[0])
       print '<td align="center">',
       if abs(n["quadrupole"])>0:
           print'[http://www.pascal-man.com/periodic-table/'+name+'.shtml '+name+']',
       print '</td>',
       print '</tr>'
       print ' '
    print '</table>'

    
