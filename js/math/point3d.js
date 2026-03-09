class Point3d {
    constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    // Copy constructor
    static copy(point) {
        return new Point3d(point.x, point.y, point.z);
    }

    // Assign operator
    assign(point) {
        this.x = point.x;
        this.y = point.y;
        this.z = point.z;
    }

    // Comparison operator
    isEqual(point) {
        return this.x === point.x && this.y === point.y && this.z === point.z;
    }

    static distance(pa, pb) {
        return Math.sqrt((pa.x - pb.x) * (pa.x - pb.x) + (pa.y - pb.y) * (pa.y - pb.y) + (pa.z - pb.z) * (pa.z - pb.z))
    }

    static translate(vp1, vp2) {
        if (vp1 instanceof Vector && vp2 instanceof Point3d)
            return new Point3d(vp1.vx + vp2.x, vp1.vy + vp2.y, vp1.vz + vp2.z);

        if (vp1 instanceof Point3d && vp2 instanceof Vector)
            return new Point3d(vp1.x + vp2.vx, vp1.y + vp2.vy, vp1.z + vp2.vz);

        throw new Error("Point3d.translate expects (Vector, Point3d) or (Point3d, Vector)");
    }
}
