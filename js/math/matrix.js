class Matrix {
    constructor(a11 = 1, a12 = 0, a13 = 0, a21 = 0, a22 = 1, a23 = 0, a31 = 0, a32 = 0, a33 = 1) {
        // Constructor with three row vectors
        if (a11 instanceof Vector && a12 instanceof Vector && a13 instanceof Vector) {
            this.a11 = a11.vx; this.a12 = a11.vy; this.a13 = a11.vz;
            this.a21 = a12.vx; this.a22 = a12.vy; this.a23 = a12.vz;
            this.a31 = a13.vx; this.a32 = a13.vy; this.a33 = a13.vz;
            return;
        }

        // Constructor with 9 numeric components
        if (typeof a11 === "number" && typeof a12 === "number" && typeof a13 === "number" &&
            typeof a21 === "number" && typeof a22 === "number" && typeof a23 === "number" &&
            typeof a31 === "number" && typeof a32 === "number" && typeof a33 === "number") {
            this.a11 = a11; this.a12 = a12; this.a13 = a13;
            this.a21 = a21; this.a22 = a22; this.a23 = a23;
            this.a31 = a31; this.a32 = a32; this.a33 = a33;
            return;
        }

        throw new Error("Matrix constructor expects either (Vector, Vector, Vector) or 9 numeric components");
    }

    // Copy constructor
    static copy(matrix) {
        return new Matrix(matrix.a11, matrix.a12, matrix.a13,
            matrix.a21, matrix.a22, matrix.a23,
            matrix.a31, matrix.a32, matrix.a33);
    }

    static multiply(mn1, mn2) {
        if (mn1 instanceof Matrix && mn2 instanceof Matrix) {
            return new Matrix(
                mn1.a11 * mn2.a11 + mn1.a12 * mn2.a21 + mn1.a13 * mn2.a31,
                mn1.a11 * mn2.a12 + mn1.a12 * mn2.a22 + mn1.a13 * mn2.a32,
                mn1.a11 * mn2.a13 + mn1.a12 * mn2.a23 + mn1.a13 * mn2.a33,
                mn1.a21 * mn2.a11 + mn1.a22 * mn2.a21 + mn1.a23 * mn2.a31,
                mn1.a21 * mn2.a12 + mn1.a22 * mn2.a22 + mn1.a23 * mn2.a32,
                mn1.a21 * mn2.a13 + mn1.a22 * mn2.a23 + mn1.a23 * mn2.a33,
                mn1.a31 * mn2.a11 + mn1.a32 * mn2.a21 + mn1.a33 * mn2.a31,
                mn1.a31 * mn2.a12 + mn1.a32 * mn2.a22 + mn1.a33 * mn2.a32,
                mn1.a31 * mn2.a13 + mn1.a32 * mn2.a23 + mn1.a33 * mn2.a33
            );
        }

        if (typeof mn1 === "number" && mn2 instanceof Matrix) {
            return new Matrix(
                mn1 * mn2.a11, mn1 * mn2.a12, mn1 * mn2.a13,
                mn1 * mn2.a21, mn1 * mn2.a22, mn1 * mn2.a23,
                mn1 * mn2.a31, mn1 * mn2.a32, mn1 * mn2.a33);
        }

        if (mn1 instanceof Matrix && typeof mn2 === "number") {
            return new Matrix(
                mn2 * mn1.a11, mn2 * mn1.a12, mn2 * mn1.a13,
                mn2 * mn1.a21, mn2 * mn1.a22, mn2 * mn1.a23,
                mn2 * mn1.a31, mn2 * mn1.a32, mn2 * mn1.a33);
        }

        throw new Error("Matrix.multiply expects (Matrix, Matrix), (number, Matrix), or (Matrix, number)");
    }

    static multiply_by_number(matrix, number) {
        return new Matrix(
            number * matrix.a11, number * matrix.a12, number * matrix.a13,
            number * matrix.a21, number * matrix.a22, number * matrix.a23,
            number * matrix.a31, number * matrix.a32, number * matrix.a33
        );
    }

    static inverse(matrix) {
        const det = matrix.a11 * (matrix.a22 * matrix.a33 - matrix.a32 * matrix.a23) -
            matrix.a12 * (matrix.a21 * matrix.a33 - matrix.a31 * matrix.a23) +
            matrix.a13 * (matrix.a21 * matrix.a32 - matrix.a31 * matrix.a22);

        if (relativeDifference(det, 0) < REL_TOLERANCE) {
            throw new Error("Matrix is not invertible");
        }

        const M = new Matrix(
            matrix.a22 * matrix.a33 - matrix.a32 * matrix.a23,
            -(matrix.a12 * matrix.a33 - matrix.a32 * matrix.a13),
            matrix.a12 * matrix.a23 - matrix.a22 * matrix.a13,
            -(matrix.a21 * matrix.a33 - matrix.a31 * matrix.a23),
            matrix.a11 * matrix.a33 - matrix.a31 * matrix.a13,
            -(matrix.a11 * matrix.a23 - matrix.a21 * matrix.a13),
            matrix.a21 * matrix.a32 - matrix.a31 * matrix.a22,
            -(matrix.a11 * matrix.a32 - matrix.a31 * matrix.a12),
            matrix.a11 * matrix.a22 - matrix.a21 * matrix.a12
        );

        return Matrix.multiply_by_number(M, 1.0 / det);
    }

    static multiply_by_point(matrix, point) {
        return new Point3d(matrix.a11 * point.x + matrix.a12 * point.y + matrix.a13 * point.z,
            matrix.a21 * point.x + matrix.a22 * point.y + matrix.a23 * point.z,
            matrix.a31 * point.x + matrix.a32 * point.y + matrix.a33 * point.z);
    }
}
