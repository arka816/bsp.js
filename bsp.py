#!/usr/bin/env python

'''
    implements BSP tree and basic set operations using merge operations
    based on "Merging BSP Trees Yields Polyhedral Set Operations" by Naylor et. al. circa 1990

    we are not using numpy for vectorization since this code will be anyways converted to C/C++ for use with WebAssembly
'''

from __future__ import annotations # PEP 563: postponed evaluation of annotations

from algebra import Vector3d

__author__ = 'arka'
__copyright__ = ""

__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Arkaprava Ghosh"
__email__ = "arkaprava.mail@gmail.com"
__status__ = "Development"

# Universe Box Dimensions (in millimetres ~ 1 km)
# preferably 2 * the largest representable floating point value to ensure no arithmetic underflow/overflow
DIM_UNIVERSE = 1e7


class HyperPlane3D:
    '''
        implements a hyperplane embedded in 3 dimensions

    '''
    def __init__(self, normal : Vector3d, intercept : int) -> None:
        self.normal = normal
        self.intercept = intercept

    def _check_side(self, point : Vector3d) -> int:
        '''
            checks if a point is outside or inside a hyperplane
            by default: if ap_x + bp_y + cp_z > d : in positive half-space (1)
                        else: in negative half-space(0)
        '''
        return int(self.normal.dot(point) > self.intercept)

    def _intercept(self, point : Vector3d) -> tuple(float, bool):
        '''
            to introduce numerical robustness, we consider three cases the point can be w.r.t. the hyperplane
            - more than ε distance from the plane on its negative side
            - more than ε distance from the plane on its positive side
            - less than ε distance from the plane on any side
        '''
        intercept = self.normal.dot(point) 
        intercept_code = 1 if intercept > self.intercept else (-1 if intercept < self.intercept else 0)
        return (intercept, intercept_code)



class BSPNode:
    def __init__(self, leaf : bool, isPosHS : bool, hp : HyperPlane3D | None = None) -> None:
        '''
            leaf: if node is a leaf; by default left child leaf nodes are OUT nodes and right child lead nodes are IN nodes
            hp  : binary partitioner hyperplane (unbounded)
        '''
        self.leaf = leaf
        self.bp = {
            'hp': hp,
            'hp_brep': self._convert_hp_to_brep()
        }

        self.parent = None
        self.neg = None
        self.pos = None
        self.isPosHS = None


    def _convert_hp_to_brep(self) -> list[Vector3d]:
        '''
            converts hyperplane of binary partitioner into a brep representation
            since B-reps cannot represent unbounded spaces, we clip the hyperplane
            at the peripheri of the Universe Box
        '''
        # find the universe box surface with the closest normal
        # normals being : [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]
        normal = self.bp['hp'].normal
        planeIndex = max(range(3), key=normal.__getitem__)

        a, b = [normal[i] for i in range(3) if i != planeIndex ]
        c = normal[planeIndex]
        d = self.bp['hp'].intercept
        
        # find the projection of this plane on the hyperplane of the bp to get the brep for hp
        # orthogonally project along the x_i axis the x_i = 0 plane onto the bp
        # plane: ax + by + cz = d
        # while projecting any arbitrary plane k, solve for k:
        # for i, j = M, M; M, -M; -M, M; -M, -M
        brep = []
        for x_i in [-DIM_UNIVERSE, DIM_UNIVERSE]:
            for x_j in [-DIM_UNIVERSE, DIM_UNIVERSE]:
                x_k = (d - (a * x_i + b * x_j)) / c
                brep.append(Vector3d(x_i, x_j, x_k))

        return brep


    def _gen_shp(self) -> list[Vector3d]:
        '''
            converts brep of hyperplane of bp into a polygonal sub-hyperplane
            such that bp.shp = bp.hp ∩ bp.region (or the region that the bp of a particular node partitions)
            passes the b-rep of the hp up through the tree upto the root to generate shp

            if bp.hp_brep does not intersect the hp of any of it's parent nodes 
            => bp.hp does not partition bp.region at all
            => bsp tree is probably malformed

            to check if an edge of a polygon is intersected by a plane:
                check if its two end vertices are on two opposite sides of the hyperplane
        '''

        brep = self.bp['hp_brep']
        node = self
        while node is not None and node.parent is not None:
            brep = partition_brep_by_bp(brep, node.parent.bp['hp'], node.isPosHS)

        self.bp['shp_brep'] = brep

'''
    ε - algebra subroutines
    - check parallel : 
        - two defintions :
            - if two hyperplanes (or sub-hyperplanes) do not intersect within universe box, they are parallel
                this uses determinants (if det value is less than epsilon, then it is effectively zero)
            - sine of the angle between the two normals is less than a particular epsilon
'''

def check_parallel(bp1 : dict, bp2 : dict, defn : int, shp : bool) -> bool:
    '''
        defn : 1 / 2
        shp  : use the entire hyperplane or the sub-hyperplane polygon
    '''
    if defn == 1:
        if shp:
            pass
        else:
            n_1, n_2 = bp1.hp.normal, bp2.hp.normal





def partition_brep_by_bp(brep : list[Vector3d], hp : HyperPlane3D, isPosHS: int, geom : bool = False) -> list[Vector3d] | int:
    '''
        Bi-Partition_Bps

        determine the geometric configuration of the shp/hp b-rep polygon and the partitioning hyperplane
        i.e. determine the intersection P1.shp ∩ P2.hp

        partition a b-rep polygon using a hyperplane
        and only return the polygon inside the specified half-space

        isPosHs: whether to consider the half of the brep in the positive halfspace or negative of the hp

        geom: whether to just determine geometry or partition the b-rep
    '''
    intercepts = [hp._intercept(point) for point in brep]
    intercepts.append(intercepts[0])

    if geom:
        # just determine geometric configuration (1 of 7 configs)
        pass
    else:
        # if the intercepts of the two vertices of an edge are on two different sides of the hyperplane intercept
        # the hyperplane intersects the edge
        i = 0
        while i < len(intercepts) - 1:
            u, v = intercepts[i], intercepts[i+1]
            if u[1] != v[1]:
                if u[1] == 0 or v[1] == 0:
                    # one of the vertices lies on the hyperplane (both ca)
                    pass
                else:
                    # edge (i -> i+1) is intersected
                    intersection = (brep[i] * abs(u[0] - hp.intercept) + brep[i+1] * abs(v[0] - hp.intercept)) \
                        / (abs(u[0] - hp.intercept) + abs(v[0] - hp.intercept))
                    intercepts.insert(i + 1, (hp.intercept, isPosHS))
                    brep.insert(i + 1, intersection)
                    i += 1
            i += 1

        brep_hs = [brep[i] for i in range(len(brep)) if intercepts[i][1] == isPosHS]
        return brep_hs


# def partition_brep_by_brep(p1 : list[Vector3d], p2 : list[Vector3d]) -> 

def mutual_geometry_polygon(p1 : list[Vector3d], p2 : list[Vector3d]) -> int:
    '''
        determines the mutual geometric configuration of two sub-hyperplanes
        represented as b-rep polygons
    '''
    pass


def partition_bspt_by_bp(bspTree : BSPNode, bp : dict) -> tuple(BSPNode, BSPNode):
    '''
        implementation of Partition_Bspt

        partition a bsp tree T by a binary partitioner (actually the shp of a bp)
        into T- and T+ such that:
        T- = T ∩ P.hs-
        T+ = T ∩ P.hs+
    '''
    # TODO: might need to clip the sub-hyperplane of the binary partitioner to the domain of bspTree.root
    # probably won't; since bspTree.root would be partitioned by bp in the previous step

    # step 1: check the mutual geometrical position of bspTree.bp.shp and bp.shp - 7 cases
    return None, None
    

def invert_tree(T : BSPNode) -> BSPNode:
    '''
        invert tree : convert IN nodes to OUT nodes and vice versa
        since left nodes are by default OUT and right nodes by default OUT we need to exchange the child nodes for each node
    '''
    if T.leaf:
        return T.leaf

    neg_subtree = invert_tree(T.pos)
    pos_subtree = invert_tree(T.neg)
    T.neg, T.pos = neg_subtree, pos_subtree
    return T


def tree_op_cell(T1 : BSPNode, T2 : BSPNode, op : str) -> BSPNode:
    '''
        performs set operation on a tree and a cell or on two cells as defined by op argument
    '''
    if not T2.leaf:
        # make T2 the leaf for convenience
        T1, T2 = T2, T1

    if op == '∩':
        if T1.isPosHS:
            # IN node
            return BSPNode(True, True)
        else:
            # OUT node
            return T2
    elif op == '∪':
        if T1.isPosHS:
            # IN node
            return T2
        else:
            # OUT node
            return BSPNode(True, False)
    elif op == '-':
        if T1.isPosHS:
            # IN node
            return invert_tree(T2)
        else:
            # OUT node
            return BSPNode(True, False)


def merge_trees(T1 : BSPNode, T2 : BSPNode, op : str) -> BSPNode:
    '''
        implementation of Merge_Bspts

        merges two BSP trees recursively into one single tree based on the set operation specified
    '''
    if T1.leaf or T2.leaf:
        return tree_op_cell(T1, T2, op)
    else:
        # T2_HS_neg = T2 ∩ T1.bp.hs-
        # T2_HS_pos = T2 ∩ T1.bp.hs+
        T2_HS_neg, T2_HS_pos = partition_bspt_by_bp(T2, T1.bp)
        
        val = BSPNode(False)
        neg_subtree = merge_trees(T1.neg, T2_HS_neg)
        pos_subtree = merge_trees(T1.pos, T2_HS_pos)

        neg_subtree.parent = val
        pos_subtree.parent = val

        neg_subtree.isPosHS = False
        pos_subtree.isPosHS = True

        val.neg = neg_subtree
        val.pos = pos_subtree

        val.bp = T1.bp

        return val
