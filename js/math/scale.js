class Scale {
    // Constructor for surface scaling.
    // (xe0, ye0) - top left corner of the screen area
    // ewidth, eheight - width and height of the screen area
    // (xr0, yr0) - center of the real area,
    // rwidth, wheight - width and height of the real area
    constructor(xe0, ye0, ewidth, eheight, xr0, yr0, rwidth, wheight) {
        this.A = ewidth / rwidth;
        this.B = xe0 - this.A * (xr0 - rwidth / 2);
        this.C = -eheight / wheight;
        this.D = ye0 - this.C * (yr0 + wheight / 2);

        this.E = rwidth / ewidth;
        this.F = xr0 - rwidth / 2 - this.E * xe0;
        this.G = -wheight / eheight;
        this.H = yr0 + wheight / 2 - this.G * ye0;
    }

    get_screen_x(x) {
        return Math.floor(this.A * x + this.B);
    }

    get_screen_y(y) {
        return Math.floor(this.C * y + this.D);
    }

    get_real_x(xe) {
        return this.E * xe + this.F;
    }

    get_real_y(ye) {
        return this.G * ye + this.H;
    }
}
