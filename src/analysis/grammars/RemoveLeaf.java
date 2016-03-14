package vdaoengine.analysis.grammars;

import java.util.Vector;
import vdaoengine.data.*;

public class RemoveLeaf extends GrammarOperation {

	public boolean applyToElement(Graph graph, VDataSet data, Vector taxons, GraphElement el) {
		boolean applicable = false;
		if(el instanceof Node){
			Node n = (Node)el;
			//if(graph.starLeaves.get(n.key)!=null)if(((Vector)graph.outgoingEdges.get(n.key)).size()==1){
			if(((Vector)graph.outgoingEdges.get(n.key)).size()==1){
				applicable = true;
				Vector v = (Vector)graph.starLeaves.get(n.key);
				Star s = (Star)v.get(0);
				int k=-1;
				for(int i=0;i<s.neighbours.size();i++)
					if(((Node)s.neighbours.get(i)).key.equals(n.key))
						k = i;
				s.neighbours.remove(k);
				s.elasticity = graph.StarDefaultElasticity[s.neighbours.size()];
				Vector oe = (Vector)graph.outgoingEdges.get(n.key);
				for(int i=0;i<oe.size();i++) graph.removeEdge((Edge)oe.get(i));
				graph.removeNode(n.key);
				
				// If the star became an edge
				if(s.neighbours.size()==1){
					//graph.addEdge(s.centralNode, (Node)s.neighbours.get(0));
					graph.removeStar(s);
				}
				
			}
		}
		return applicable;
	}

}
