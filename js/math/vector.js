class Vector {
    constructor(vx = 0, vy = 0, vz = 0) {
        if (vx instanceof Point3d && vy instanceof Point3d) {
            // Constructor with two points
            const begPoint = vx;
            const endPoint = vy;
            this.vx = endPoint.x - begPoint.x;
            this.vy = endPoint.y - begPoint.y;
            this.vz = endPoint.z - begPoint.z;
            return;
        }

        // Constructor with components
        if (typeof vx === "number" && typeof vy === "number" && typeof vz === "number") {
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            return;
        }

        throw new Error("Vector constructor expects either (Point3d, Point3d) or (number, number, number)");
    }

    // Copy constructor
    static copy(vector) {
        return new Vector(vector.vx, vector.vy, vector.vz);
    }

    // Assign operator
    assign(vector) {
        this.vx = vector.vx;
        this.vy = vector.vy;
        this.vz = vector.vz;
        return this;
    }

    // Comparison operator
    isEqual(vector) {
        return this.vx === vector.vx && this.vy === vector.vy && this.vz === vector.vz;
    }

    // Negation operator
    negate() {
        this.vx = -this.vx;
        this.vy = -this.vy;
        this.vz = -this.vz;
        return this;
    }

    // Length function
    length() {
        return Math.sqrt(this.vx * this.vx + this.vy * this.vy + this.vz * this.vz);
    }

    normalize() {
        const len = this.length();
        if (len > 0) {
            this.vx /= len;
            this.vy /= len;
            this.vz /= len;
        }
        return this;
    }

    // Adding two vectors
    static add(vector1, vector2) {
        return new Vector(vector1.vx + vector2.vx, vector1.vy + vector2.vy, vector1.vz + vector2.vz);
    }

    // Product of a vector by a number
    multiplyByNumber(number) {
        return new Vector(this.vx * number, this.vy * number, this.vz * number);
    }

    // Dot product of two vectors
    static dot_product(vector1, vector2) {
        return vector1.vx * vector2.vx + vector1.vy * vector2.vy + vector1.vz * vector2.vz;
    }

    // Cross product of two vectors
    static cross_product(vector1, vector2) {
        return new Vector(
            vector1.vy * vector2.vz - vector1.vz * vector2.vy,
            vector1.vz * vector2.vx - vector1.vx * vector2.vz,
            vector1.vx * vector2.vy - vector1.vy * vector2.vx
        );
    }
}
