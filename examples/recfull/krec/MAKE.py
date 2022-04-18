## Compile ##
import makefile as mk

mk.myl['utils'] = True
mk.myl['nldd'] = True
mk.myl['anafull'] = True

mk.makefile(myl=mk.myl)
mk.compile(remove=False)

