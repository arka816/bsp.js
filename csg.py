class CSGNode:
    '''
        implementation of half-space CSG tree
        leaves: non-empty regular half-spaces
        nodes: regularized set operations on said half-space leaves or other nodes
    '''
    def __init__(self, leaf, op=None, surface=None):
        self.leaf = leaf
        
        if self.leaf:
            assert(surface is not None)
            self.surface = surface
        else:
            assert()
            self.op = op