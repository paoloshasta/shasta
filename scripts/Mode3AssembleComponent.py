#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Load a mode3::AssemblyGraph representing a connected component of the primary graph and assemble it.')
parser.add_argument('component', type=int, 
    help='The connected component to assemble.')
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)  
        
arguments = parser.parse_args() 



a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()
shasta.openPerformanceLog('Mode3AssembleComponent.log')
fileName = 'AssemblyGraph-' + str(arguments.component) + '.data'
a.mode3AssembleComponent(fileName, 0, arguments.debug)
 