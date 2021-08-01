"use strict";

class HeatMethod {
	/**
	 * This class implements the {@link http://cs.cmu.edu/~kmcrane/Projects/HeatMethod/ heat method} to compute geodesic distance
	 * on a surface mesh.
	 * @constructor module:Projects.HeatMethod
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} F The mean curvature flow operator built on the input mesh.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);

		// TODO: build laplace and flow matrices
		this.A = geometry.laplaceMatrix(this.vertexIndex); // placeholder

		// mean curvature flow
		let M = geometry.massMatrix(this.vertexIndex);
		let h = this.geometry.meanEdgeLength();

		// t = h*h
		let F = M.plus(this.A.timesReal(h*h));

		this.F = F; // placeholder
	}

	/**
	 * Computes the vector field X = -∇u / |∇u|.
	 * @private
	 * @method module:Projects.HeatMethod#computeVectorField
	 * @param {module:LinearAlgebra.DenseMatrix} u A dense vector (i.e., u.nCols() == 1) representing the
	 * heat that is allowed to diffuse on the input mesh for a brief period of time.
	 * @returns {Object} A dictionary mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 */
	computeVectorField(u) {
		// TODO
		let X = {};

		for (let f of this.geometry.mesh.faces) {
			let N = this.geometry.faceNormal(f);
			let Af= this.geometry.area(f);
			let du= new Vector();
	
			for (let h of f.adjacentHalfedges()) {
				let ui = u.get(this.vertexIndex[h.prev.vertex], 0);

				// ui * (N \times ei)
				// let dui= N.cross(this.geometry.vector(h));
				let dui = this.geometry.vector(h);
				dui.scaleBy(ui);

				du.decrementBy(dui);
			}

			du = N.cross(du);
			du.divideBy(2.0*Af);

			X[f] = du.unit();
		}

		return X; // placeholder
	}

	/**
	 * Computes the integrated divergence ∇.X.
	 * @private
	 * @method module:Projects.HeatMethod#computeDivergence
	 * @param {Object} X The vector field -∇u / |∇u| represented by a dictionary
	 * mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeDivergence(X) {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		let div = DenseMatrix.zeros(vertices.length, 1);

		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let d = 0.0;
			
			for (let h of v.adjacentHalfedges()) {
				if (!h.onBoundary) {
					let e = this.geometry.vector(h);

					let Xj = X[h.face];
					d += this.geometry.cotan(h) * Xj.dot(e);
					
					h = h.twin;
					Xj = X[h.face];
					d += this.geometry.cotan(h) * Xj.dot(e);
				}
			}
			
			div.set(0.5*d, i, 0);
		}

		return div; // placeholder
	}

	/**
	 * Shifts φ such that its minimum value is zero.
	 * @private
	 * @method module:Projects.HeatMethod#subtractMinimumDistance
	 * @param {module:LinearAlgebra.DenseMatrix} phi The (minimum 0) solution to the poisson equation Δφ = ∇.X.
	 */
	subtractMinimumDistance(phi) {
		let min = Infinity;
		for (let i = 0; i < phi.nRows(); i++) {
			min = Math.min(phi.get(i, 0), min);
		}

		for (let i = 0; i < phi.nRows(); i++) {
			phi.set(phi.get(i, 0) - min, i, 0);
		}
	}

	/**
	 * Computes the geodesic distances φ using the heat method.
	 * @method module:Projects.HeatMethod#compute
	 * @param {module:LinearAlgebra.DenseMatrix} delta A dense vector (i.e., delta.nCols() == 1) containing
	 * heat sources, i.e., u0 = δ(x).
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	compute(delta) {
		// TODO
		// Choleksy factorization
		let Fllt = this.F.chol();
		let Allt = this.A.chol();

		// integrate the heat flow => u
		let u = Fllt.solvePositiveDefinite(delta);

		// Evaluate the vector field => X
		let X = this.computeVectorField(u);

		// Solve the Poisson equation => phi
		// note that our laplacian is positive semidefinite, the equation should be: -L phi = div
		let div = this.computeDivergence(X);
		let phi = Allt.solvePositiveDefinite(div.negated());

		// since φ is unique up to an additive constant, it should
		// be shifted such that the smallest distance is zero
		this.subtractMinimumDistance(phi);

		return phi;
	}
}
