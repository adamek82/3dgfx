class Point2d {
    constructor(x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }

    // Copy constructor
    static copy(point) {
        return new Point2d(point.x, point.y);
    }

    // Assign operator
    assign(point) {
        this.x = point.x;
        this.y = point.y;
    }

    // Comparison operator
    isEqual(point) {
        return this.x === point.x && this.y === point.y;
    }

    static distance(pa, pb) {
        return Math.sqrt((pa.x - pb.x) * (pa.x - pb.x) + (pa.y - pb.y) * (pa.y - pb.y))
    }
}
