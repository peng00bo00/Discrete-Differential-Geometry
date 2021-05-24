"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // TODO
                // assign vertices
                for (let i = 0; i < mesh.vertices.length; i++) {
                        mesh.vertices[i].index = i;
                }

                // assign edges
                for (let i = 0; i < mesh.edges.length; i++) {
                        mesh.edges[i].index = i;
                }

                // assign faces
                for (let i = 0; i < mesh.faces.length; i++) {
                        mesh.faces[i].index = i;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                // TODO
                let T = new Triplet(mesh.edges.length, mesh.vertices.length);

                for (const edge of mesh.edges) {
                        T.addEntry(1, edge.index, edge.halfedge.vertex.index);
                        T.addEntry(1, edge.index, edge.halfedge.twin.vertex.index);
                }

                return SparseMatrix.fromTriplet(T);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                // TODO
                let T = new Triplet(mesh.faces.length, mesh.edges.length);

                for (const face of mesh.faces) {
                        let edge1 = face.halfedge.edge;
                        let edge2 = face.halfedge.next.edge;
                        let edge3 = face.halfedge.next.next.edge;

                        T.addEntry(1, face.index, edge1.index);
                        T.addEntry(1, face.index, edge2.index);
                        T.addEntry(1, face.index, edge3.index);
                }

                return SparseMatrix.fromTriplet(T);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                // TODO
                let V = DenseMatrix.zeros(this.mesh.vertices.length, 1);

                for (const v of subset.vertices) {
                        V.set(1, v, 0);
                }

                return V;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                // TODO
                let E = DenseMatrix.zeros(this.mesh.edges.length, 1);

                for (const e of subset.edges) {
                        E.set(1, e, 0);
                }

                return E;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                // TODO
                let F = DenseMatrix.zeros(this.mesh.faces.length, 1);

                for (const f of subset.faces) {
                        F.set(1, f, 0);
                }

                return F;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
		let meshsub = MeshSubset.deepCopy(subset);

                // add adjacent edges
		let E0 = this.A0;
		let V = this.buildVertexVector(meshsub);
		let adjEdges = E0.timesDense(V);

		for (let i = 0; i < adjEdges.nRows(); i++) {
			if (adjEdges.get(i, 0) > 0) {
                                let idx = this.mesh.edges[i].index;
				meshsub.addEdge(idx);
			}
		}

                // add adjacent faces
		let E1 = this.A1;
		let E = this.buildEdgeVector(meshsub);
		let adjFaces = E1.timesDense(E);

		for (let i = 0; i < adjFaces.nRows(); i++) {
			if (adjFaces.get(i, 0) > 0) {
                                let idx = this.mesh.faces[i].index;
				meshsub.addFace(idx);
			}
		}

		return meshsub;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO
                let meshsub = MeshSubset.deepCopy(subset);

                // add adjacent edges
                let E1 = this.A1;
		let F = this.buildFaceVector(meshsub);
		let adjEdges = E1.transpose().timesDense(F);

                for (let i = 0; i < adjEdges.nRows(); i++) {
			if (adjEdges.get(i, 0) > 0) {
                                let idx = this.mesh.edges[i].index;
				meshsub.addEdge(idx);
			}
		}

                // add adjacent vertices
                let E0 = this.A0;
                let E = this.buildEdgeVector(meshsub);
                let adjVertices = E0.transpose().timesDense(E);

                for (let i = 0; i < adjVertices.nRows(); i++) {
			if (adjVertices.get(i, 0) > 0) {
                                let idx = this.mesh.vertices[i].index;
				meshsub.addVertex(idx);
			}
		}

                return meshsub; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO
                // link = closure(star) - star(closure)
                let closure_star = this.closure(this.star(subset));
		let star_closure = this.star(this.closure(subset));

		closure_star.deleteSubset(star_closure);

		return closure_star;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                if (subset.equals(this.closure(subset))) {
			return true;
		}

                return false;
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
		if (!this.isComplex(subset)) {
			return -1;
		}

		let pureComplex = new MeshSubset;
		if (subset.faces.size > 0) {
			for (let f of subset.faces) {
				pureComplex.addFace(f);
			}

			if (this.closure(pureComplex).equals(subset)) {
				return 2;
			} else {
				return -1;
			}
		} else if (subset.edges.size > 0) {
			for (let e of subset.edges) {
				pureComplex.addEdge(e);
			}

			if (this.closure(pureComplex).equals(subset)) {
				return 1;
			} else {
				return -1;
			}
		} else if (subset.vertices.size > 0) {
			return 0;
		}

		return -1;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO
		let myBoundary = new MeshSubset;

		if (subset.faces.size > 0) {
			let FEmatrix = this.A1.transpose();
			let Fvector = this.buildFaceVector(subset);
			let adjEdges = FEmatrix.timesDense(Fvector);

			for (let i = 0; i < adjEdges.nRows(); i++) {
				if (adjEdges.get(i, 0) == 1) {
					myBoundary.addEdge(i);
				}
			}

		} else if (subset.edges.size > 0) {
			let EVmatrix = this.A0.transpose();
			let Evector = this.buildEdgeVector(subset);
			let adjVertices = EVmatrix.timesDense(Evector);

			for (let i = 0; i < adjVertices.nRows(); i++) {
				if (adjVertices.get(i, 0) == 1) {
					myBoundary.addVertex(i);
				}
			}
		}

		return this.closure(myBoundary);
        }
}
