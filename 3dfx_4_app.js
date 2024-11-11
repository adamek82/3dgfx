// Define a small tolerance value for relative comparison
const REL_TOLERANCE = 1e-10;

// Function to calculate relative difference
function relativeDifference(a, b) {
    let absA = Math.abs(a);
    let absB = Math.abs(b);
    let maxAbs = Math.max(absA, absB);
    return maxAbs === 0 ? 0 : Math.abs(a - b) / maxAbs;
}

class Point2d {
    constructor (x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }

    // Copy constructor
    static copy(point) {
        return new Point3d(point.x, point.y);
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
        return Math.sqrt((pa.x - pb.x)*(pa.x - pb.x) + (pa.y - pb.y)*(pa.y - pb.y))
    }
}

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
        return Math.sqrt((pa.x - pb.x)*(pa.x - pb.x) + (pa.y - pb.y)*(pa.y - pb.y) + (pa.z - pb.z)*(pa.z - pb.z))
    }

    static translate(vp1, vp2) {
        if (vp1 instanceof Vector && vp2 instanceof Point3d)
            return new Point3d(vp1.vx + vp2.x, vp1.vy + vp2.y + vp1.vz + vp2.z);
        else if (vp1 instanceof Point3d && vp2 instanceof Vector)
            return new Point3d(vp1.x + vp2.vx, vp1.y + vp2.vy, vp1.z + vp2.vz);
        else
            return new Point3d();     // TODO: consider throwing an exception here
    }
}

class Vector {
    constructor(vx = 0, vy = 0, vz = 0) {
        if (vx instanceof Point3d && vy instanceof Point3d) {
            // Constructor with two points
            const begPoint = vx;
            const endPoint = vy;
            this.vx = endPoint.x - begPoint.x;
            this.vy = endPoint.y - begPoint.y;
            this.vz = endPoint.z - begPoint.z;
        } else {
            // Constructor with components
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
        }
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

class Matrix {
    constructor(a11 = 1, a12 = 0, a13 = 0, a21 = 0, a22 = 1, a23 = 0, a31 = 0, a32 = 0, a33 = 1) {
        if (a11 instanceof Vector && a12 instanceof Vector && a13 instanceof Vector) {
            this.a11 = a11.vx; this.a12 = a11.vy; this.a13 = a11.vz;
            this.a21 = a12.vx; this.a22 = a12.vy; this.a23 = a12.vz;
            this.a31 = a13.vx; this.a32 = a13.vy; this.a33 = a13.vz;
        }
        else {
            this.a11 = a11; this.a12 = a12; this.a13 = a13;
            this.a21 = a21; this.a22 = a22; this.a23 = a23;
            this.a31 = a31; this.a32 = a32; this.a33 = a33;
        }
    }

    // Copy constructor
    static copy(matrix) {
        return new Matrix(matrix.a11, matrix.a12, matrix.a13, 
                          matrix.a21, matrix.a22, matrix.a23, 
                          matrix.a31, matrix.a32, matrix.a33);
    }

    static multiply(mn1, mn2) {
        if (mn1 instanceof Matrix && mn2 instanceof Matrix) {
            return new Matrix(mn1.a11*mn2.a11 + mn1.a12*mn2.a21 + mn1.a13*mn2.a31,
                              mn1.a11*mn2.a12 + mn1.a12*mn2.a22 + mn1.a13*mn2.a32,
                              mn1.a11*mn2.a13 + mn1.a12*mn2.a23 + mn1.a13*mn2.a33,
                              mn1.a21*mn2.a11 + mn1.a22*mn2.a21 + mn1.a23*mn2.a31,
                              mn1.a21*mn2.a12 + mn1.a22*mn2.a22 + mn1.a23*mn2.a32,
                              mn1.a21*mn2.a13 + mn1.a22*mn2.a23 + mn1.a23*mn2.a33,
                              mn1.a31*mn2.a11 + mn1.a32*mn2.a21 + mn1.a33*mn2.a31,
                              mn1.a31*mn2.a12 + mn1.a32*mn2.a22 + mn1.a33*mn2.a32,
                              mn1.a31*mn2.a13 + mn1.a32*mn2.a23 + mn1.a33*mn2.a33
            );
        }
        else if (mn2 instanceof Matrix) {
            return new Matrix(mn1*mn2.a11, mn1*mn2.a12, mn1*mn2.a13, 
                              mn1*mn2.a21, mn1*mn2.a22, mn1*mn2.a23, 
                              mn1*mn2.a31, mn1*mn2.a32, mn1*mn2.a33);
        }
        else {
            return new Matrix(mn2*mn1.a11, mn2*mn1.a12, mn1*mn1.a13, 
                mn2*mn1.a21, mn2*mn1.a22, mn2*mn1.a23, 
                mn2*mn1.a31, mn2*mn1.a32, mn2*mn1.a33);
        }
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
        return new Point3d(matrix.a11*point.x + matrix.a12*point.y + matrix.a13*point.z,
                         matrix.a21*point.x + matrix.a22*point.y + matrix.a23*point.z,
                         matrix.a31*point.x + matrix.a32*point.y + matrix.a33*point.z);
    }
}

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
        
        if (relativeDifference(newp.z, 0) < REL_TOLERANCE) {
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

        if (relativeDifference(newp.z, 0) < REL_TOLERANCE) {
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

// global variables
var canvas_id = "mainCanvas",
	canvas = document.getElementById(canvas_id),
	ctx = canvas.getContext("2d"),
	width,
	height,
    R = 40,                     // You
    FI = 60,                    //
    TETA = 80,                  //
    SCR_DIST = 50,              //
    R_obs = 5,                  // Virtual observer
    FI_obs = 120,               // 
    TETA_obs = 60;              // 
    
 const e_dist = 0.5,
    e_width = 1, 
    e_height = 0.8;

// Cube vertices
const k = [
    new Point3d(1, 1, -1),  new Point3d(-1, 1, -1), new Point3d(-1, -1, -1), new Point3d(1, -1, -1),
    new Point3d(1, 1, 1),   new Point3d(-1, 1, 1),  new Point3d(-1, -1, 1),  new Point3d(1, -1, 1)
];

// helper functions
function _(d) {
	return document.getElementById(d)
}

function setPixel(x, y, r, g, b)
{
	if (typeof g !== "undefined" && typeof b !== "undefined")
		ctx.fillStyle = "rgb("+r+","+g+","+b+")";
	else
		ctx.fillStyle = r;
	ctx.fillRect(x, y, 1, 1);
}

function drawLine(x1, y1, x2, y2, r, g, b)
{
    if (typeof g !== "undefined" && typeof b !== "undefined")
        ctx.strokeStyle = "rgb("+r+","+g+","+b+")";
    else
        ctx.strokeStyle = r;
    
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.lineWidth = 1;
    ctx.stroke();
}

function sphericalToCartesian(r, phi, theta)
{
    const k = Math.PI / 180.0;      // Conversion factor from degrees to radians

    const x = r * Math.sin(k * theta) * Math.cos(k * phi);
    const y = r * Math.sin(k * theta) * Math.sin(k * phi);
    const z = r * Math.cos(k * theta);

    return new Point3d(x, y, z);
}

function get_points_3d(d, points)
{
    let points2d = [];
    for (let p of points) {
        let point2d = d.point_3d(p);
        if (!point2d)
            return null;
        points2d.push(point2d);
    }
    return points2d;
}

function draw_main_axes(d) {
    const a = R / 5;
    const px1 = new Point3d(a, 0, 0), py1 = new Point3d(0, a, 0), pz1 = new Point3d(0, 0, a);
    const px2 = new Point3d(-a, 0, 0), py2 = new Point3d(0, -a, 0), pz2 = new Point3d(0, 0, -a);

    const points = [px1, px2, py1, py2, pz1, pz2];
    const screenPoints = get_points_3d(d, points);

    if (screenPoints) {
        const axes_color = "rgb(192, 192, 192)";
        // Draw X axis
        drawLine(screenPoints[0].x, screenPoints[0].y, screenPoints[1].x, screenPoints[1].y, axes_color);
        // Draw Y axis
        drawLine(screenPoints[2].x, screenPoints[2].y, screenPoints[3].x, screenPoints[3].y, axes_color);
        // Draw Z axis
        drawLine(screenPoints[4].x, screenPoints[4].y, screenPoints[5].x, screenPoints[5].y, axes_color);

        const dr = 1, ddr = 0.1;

        // Helper function to draw scales on an axis
        function drawAxisScale(coordinate) {
            for (let i = 0; i <= a; i += dr) {
                let pointsToDraw = [
                    { p1: new Point3d(...(coordinate === 'x' ? [i, 0, -ddr] : coordinate === 'y' ? [0, i, -ddr] : [-ddr, 0, i])),
                      p2: new Point3d(...(coordinate === 'x' ? [i, 0, ddr] : coordinate === 'y' ? [0, i, ddr] : [ddr, 0, i])) },
                    { p1: new Point3d(...(coordinate === 'x' ? [-i, 0, -ddr] : coordinate === 'y' ? [0, -i, -ddr] : [-ddr, 0, -i])),
                      p2: new Point3d(...(coordinate === 'x' ? [-i, 0, ddr] : coordinate === 'y' ? [0, -i, ddr] : [ddr, 0, -i])) }
                ];

                pointsToDraw.forEach(({ p1, p2 }) => {
                    let sp1 = d.point_3d(p1);
                    let sp2 = d.point_3d(p2);

                    if (sp1 && sp2) {
                        drawLine(sp1.x, sp1.y, sp2.x, sp2.y, axes_color);
                    }
                });
            }
        }

        // Draw scales on X, Y, and Z axes
        drawAxisScale('x');
        drawAxisScale('y');
        drawAxisScale('z');
    }
}

// Drawing the coordinate system of the virtual observer
// Object d is used for drawing, d_obs is used for calculations of the visible scene's structure
function draw_observer_axes(d, d_obs, observer) {
    const center = new Point3d(0, 0, 0);
    const M = d_obs.M;  // The matrix contains the description of the observer's base
    const vx = new Vector(M.a11, M.a12, M.a13);  // Base vectors of the observer's system
    const vy = new Vector(M.a21, M.a22, M.a23);
    const vz = new Vector(M.a31, M.a32, M.a33);
    
    let px = Point3d.translate(observer, vx.multiplyByNumber(3));  // Endpoints of the observer's base vectors
    let py = Point3d.translate(observer, vy.multiplyByNumber(3));  // Aesthetically extended 3 times
    let pz = Point3d.translate(observer, vz.multiplyByNumber(3));

    let xo, yo, x0, y0, xc, yc, x1, y1, x2, y2, x3, y3;  // Screen equivalents

    // Color for the observer's system
    const observer_system_color = "rgb(192, 192, 192)";

    if (d.point_3d(observer) && d.point_3d(new Point3d(0, 0, 0)) && 
        d.point_3d(center) && d.point_3d(px) && 
        d.point_3d(py) && d.point_3d(pz)) {
        
        let screenPoints = [observer, new Point3d(0, 0, 0), center, px, py, pz].map(p => d.point_3d(p));

        if (screenPoints.every(sp => sp)) {
            [xo, yo] = [screenPoints[0].x, screenPoints[0].y];
            [x0, y0] = [screenPoints[1].x, screenPoints[1].y];
            [xc, yc] = [screenPoints[2].x, screenPoints[2].y];
            [x1, y1] = [screenPoints[3].x, screenPoints[3].y];
            [x2, y2] = [screenPoints[4].x, screenPoints[4].y];
            [x3, y3] = [screenPoints[5].x, screenPoints[5].y];

            // Color for the observer's optical axis
            const observer_optical_axis_color = "rgb(255, 0, 0)";
            drawLine(xo, yo, xc, yc, observer_optical_axis_color);

            drawLine(xo, yo, x1, y1, observer_system_color);  // X' axis
            drawLine(xo, yo, x2, y2, observer_system_color);  // Y' axis
            drawLine(xo, yo, x3, y3, observer_system_color);  // Z' axis
        }
    }

    // Projection onto the XY plane
    pz = Point3d.copy(observer);
    pz.z = 0;
    px = Point3d.copy(pz);
    px.y = 0;
    py = Point3d.copy(pz);
    py.x = 0;

    let sp_pz = d.point_3d(pz);
    let sp_px = d.point_3d(px);
    let sp_py = d.point_3d(py);

    if (sp_pz && sp_px && sp_py) {
        [x1, y1] = [sp_pz.x, sp_pz.y];
        [x2, y2] = [sp_px.x, sp_px.y];
        [x3, y3] = [sp_py.x, sp_py.y];

        drawLine(xo, yo, x1, y1, observer_system_color);  // Line from the center to the observer
        drawLine(x0, y0, x1, y1, observer_system_color);
        drawLine(x1, y1, x2, y2, observer_system_color);
        drawLine(x1, y1, x3, y3, observer_system_color);
    }
}

//  Drawing the virtual observer's screen
//  Then, drawing the image on that small screen
//  First, the observer's projector (d_obs) is used to find metric coordinates
//  on the virtual screen.
//  Then, those coordinates are transformed by your projector (d) onto your screen
function draw_observer_screen(d, d_obs) {
    let x = new Array(8), y = new Array(8); // Pixel coordinates of the cube's vertices
    let p_screen = new Array(8);               // Full 3D points representing the screen corners
    let p;                                  // A point in your coordinate system
    let xr, yr;
    let pe = new Array(4);                  // Array to store the 2D points for the screen border

    // First, let's draw the observer's screen border
    p_screen[0] = new Point3d(-e_width / 2, e_height / 2, e_dist); // Describe the screen corners in the observer's coordinate system
    p_screen[1] = new Point3d(e_width / 2, e_height / 2, e_dist);
    p_screen[2] = new Point3d(e_width / 2, -e_height / 2, e_dist);
    p_screen[3] = new Point3d(-e_width / 2, -e_height / 2, e_dist);

    for (let i = 0; i < 4; i++) {
        p = d_obs.inverse_transform(p_screen[i]); // Transform the point back to the global coordinate system
        let screenPoint = d.point_3d(p);
        if (!screenPoint) return;            // Now present this point on your scene

        pe[i] = new Point2d(screenPoint.x, screenPoint.y);
    }

    const observer_screen_border_color = "rgb(0, 255, 0)";
    for (let i = 0; i < 4; i++) {
        drawLine(pe[i].x, pe[i].y, pe[(i + 1) % 4].x, pe[(i + 1) % 4].y, observer_screen_border_color);
    }

    for (let i = 0; i < 8; i++) {
        let screenCoords = d_obs.point_3d_real(k[i]); // Helper function gives metric coordinates
        if (!screenCoords) return;

        p = new Point3d(screenCoords.x, screenCoords.y, e_dist); // Pixel location on the virtual screen
        p = d_obs.inverse_transform(p); // Transform back to the global coordinate system

        let screenPoint = d.point_3d(p);
        if (!screenPoint) return;       // Now present this point on your scene

        x[i] = screenPoint.x;
        y[i] = screenPoint.y;
    }

    const object_image_color = "rgb(0, 0, 255)";
    drawLine(x[0], y[0], x[1], y[1], object_image_color);  // Base of the cube/pyramid
    drawLine(x[1], y[1], x[2], y[2], object_image_color);
    drawLine(x[2], y[2], x[3], y[3], object_image_color);
    drawLine(x[3], y[3], x[0], y[0], object_image_color);

    drawLine(x[4], y[4], x[5], y[5], object_image_color);  // Base of the cube
    drawLine(x[5], y[5], x[6], y[6], object_image_color);
    drawLine(x[6], y[6], x[7], y[7], object_image_color);
    drawLine(x[7], y[7], x[4], y[4], object_image_color);

    drawLine(x[0], y[0], x[4], y[4], object_image_color);  // Side edges
    drawLine(x[1], y[1], x[5], y[5], object_image_color);
    drawLine(x[2], y[2], x[6], y[6], object_image_color);
    drawLine(x[3], y[3], x[7], y[7], object_image_color);
}

function draw_cube(d)
{
    // Transform points (3d -> 2d)
    let points2d = get_points_3d(d, k);

    if (!points2d)
        return;

    // Set the color to yellow for the lines
    const yellow = "rgb(255, 255, 0)";

    // Draw the base of the cube
    drawLine(points2d[0].x, points2d[0].y, points2d[1].x, points2d[1].y, yellow);
    drawLine(points2d[1].x, points2d[1].y, points2d[2].x, points2d[2].y, yellow);
    drawLine(points2d[2].x, points2d[2].y, points2d[3].x, points2d[3].y, yellow);
    drawLine(points2d[3].x, points2d[3].y, points2d[0].x, points2d[0].y, yellow);

    // Draw the top of the cube
    drawLine(points2d[4].x, points2d[4].y, points2d[5].x, points2d[5].y, yellow);
    drawLine(points2d[5].x, points2d[5].y, points2d[6].x, points2d[6].y, yellow);
    drawLine(points2d[6].x, points2d[6].y, points2d[7].x, points2d[7].y, yellow);
    drawLine(points2d[7].x, points2d[7].y, points2d[4].x, points2d[4].y, yellow);

    // Draw the vertical edges of the cube
    drawLine(points2d[0].x, points2d[0].y, points2d[4].x, points2d[4].y, yellow);
    drawLine(points2d[1].x, points2d[1].y, points2d[5].x, points2d[5].y, yellow);
    drawLine(points2d[2].x, points2d[2].y, points2d[6].x, points2d[6].y, yellow);
    drawLine(points2d[3].x, points2d[3].y, points2d[7].x, points2d[7].y, yellow);
}

function size_change(d)
{
	_("size_info").innerHTML = d;
	w = d + "px";
	_(canvas_id).style.width = w;
	_(canvas_id).style.height = 3 * d / 4 + "px";
}

function set_zoom()
{
	width = canvas.width = canvas.clientWidth;
	height = canvas.height = canvas.clientHeight;
}

function get_settings() {
    width = canvas.width;
    height = canvas.height;

    R = _("R").value;
    _("R_info").innerHTML = R;
    FI = _("FI").value;
    _("FI_info").innerHTML = FI;
    TETA = _("TETA").value;
    _("TETA_info").innerHTML = TETA;
    SCR_DIST = _("SCR_DIST").value;
    _("SCR_DIST_info").innerHTML = SCR_DIST;

    R_obs = _("R_obs").value;
    _("R_obs_info").innerHTML = R_obs;
    FI_obs = _("FI_obs").value;
    _("FI_obs_info").innerHTML = FI_obs;
    TETA_obs = _("TETA_obs").value;
    _("TETA_obs_info").innerHTML = TETA_obs;

    draw_scene();
}

function draw_scene()
{
	set_zoom();
	
	ctx.fillStyle = 'black';
	ctx.fillRect(0, 0, width, height);

    const scr_dist_meters = SCR_DIST / 100;
    const you = sphericalToCartesian(R, FI, TETA);
    const observer = sphericalToCartesian(R_obs, FI_obs, TETA_obs);
    const pixel_width = 800, pixel_height = 600;

    let d = new C3DView(you, 0, 0, width, height, scr_dist_meters);
    let d_obs = new C3DView(observer, 0, 0, pixel_width, pixel_height, e_dist, e_width, e_height);

    draw_main_axes(d);                        // Draw the main coordinate system
    draw_observer_axes(d, d_obs, observer);   // Draw the observer's coordinate system
    draw_observer_screen(d, d_obs);           // Draw the observer's screen 
	draw_cube(d);                             // Draw the cube
}
