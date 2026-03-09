class C3DView {
    constructor(observer,                               // observer position
        xe0, ye0, ewidth, eheight,                      // screen window in pixels
        screen_dist = 0.5,                              // typical eye distance from the screen
        screen_width = 0, screen_height = 0) {          // physical description of the screen in meters

        let center = new Point3d(0, 0, 0);
        let vertical = new Vector(0, 0, 1);

        if (screen_width <= 0)                          // when no screen sizes are specified
            screen_width = ewidth * 0.0254 / 96.0;      // 0.0254 m = 1 inch. The screen has 96 pixels/inch
        if (screen_height <= 0)                         //
            screen_height = eheight * 0.0254 / 96.0;    //

        this.s = new Scale(xe0, ye0, ewidth, eheight, 0, 0, screen_width, screen_height);

        // local z' vector of the observer (indicating the center of the system)
        let v3 = new Vector(observer, center);
        v3.normalize();
        // local x' axis vector of the observer
        let v1 = Vector.cross_product(vertical, v3);
        v1.normalize();
        // ;ocal vector of the observer's y' axis
        let v2 = Vector.cross_product(v3, v1);
        // left-handedness of the observer's system
        v1.negate();
        // not necessary when the components are normalized and perpendicular
        //v2.normalize();
        // transformation matrix
        this.M = new Matrix(v1, v2, v3);

        this.f_screen_dist = screen_dist;               // parameters needed for other functions
        this.f_R = new Vector(observer, center);        // distance between the coordinate systems
    }

    point_3d(p) {
        let newp = Point3d.translate(p, this.f_R);
        newp = Matrix.multiply_by_point(this.M, newp);

        if (Math.abs(newp.z) < REL_TOLERANCE) {
            // If newp.z is too close to zero, return false to indicate invalid transformation
            return false;
        } else {
            let xe = this.s.get_screen_x(newp.x * this.f_screen_dist / newp.z);
            let ye = this.s.get_screen_y(newp.y * this.f_screen_dist / newp.z);
            return new Point2d(xe, ye);
        }
    }

    point_3d_real(p) {
        let newp = Point3d.translate(p, this.f_R);
        newp = Matrix.multiply_by_point(this.M, newp);

        if (Math.abs(newp.z) < REL_TOLERANCE) {
            return false;
        }

        let xr = newp.x * this.f_screen_dist / newp.z;
        let yr = newp.y * this.f_screen_dist / newp.z;

        return new Point2d(xr, yr);
    }

    inverse_transform(p) {
        let M_inverse = Matrix.inverse(this.M);  // Calculate the inverse of the transformation matrix
        let fR_inverse = Vector.copy(this.f_R).negate(); // Negate the translation vector
        let transformed_point = Matrix.multiply_by_point(M_inverse, p);  // Apply the inverse matrix to the point
        transformed_point = Point3d.translate(transformed_point, fR_inverse);  // Add the negated translation vector
        return transformed_point;  // Return the transformed point
    }
}
