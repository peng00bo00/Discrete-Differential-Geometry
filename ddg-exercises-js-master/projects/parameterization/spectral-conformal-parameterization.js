"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {
		// TODO
		let V = this.geometry.mesh.vertices.length;

		// ED
		let ED = this.geometry.complexLaplaceMatrix(this.vertexIndex);
		ED = ED.timesComplex(new Complex(0.5, 0.));

		// A
		let T = new ComplexTriplet(V, V);
		for (let f of this.geometry.mesh.boundaries) {
			for (let h of f.adjacentHalfedges()) {
				if (h.onBoundary) {
					let i = this.vertexIndex[h.vertex];
					let j = this.vertexIndex[h.twin.vertex];
	
					T.addEntry(new Complex(0., 0.25), i, j);
					T.addEntry(new Complex(0.,-0.25), j, i);
				}
			}
		}

		let A = ComplexSparseMatrix.fromTriplet(T);

		return ED.minus(A);
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		let flattening = this.geometry.positions; // placeholder

		let x = Solvers.solveInversePowerMethod(this.buildConformalEnergy());

		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let p = flattening[v];

			p.x = x.get(i, 0).re;
			p.y = x.get(i, 0).im;
			p.z = 0.;
		}

		// normalize flattening
		normalize(flattening, vertices);

		return flattening;
	}
}
