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

ctx.lineWidth = 1;
ctx.lineCap = "butt";
ctx.lineJoin = "miter";

function get_points_3d(d, points) {
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
                    {
                        p1: new Point3d(...(coordinate === 'x' ? [i, 0, -ddr] : coordinate === 'y' ? [0, i, -ddr] : [-ddr, 0, i])),
                        p2: new Point3d(...(coordinate === 'x' ? [i, 0, ddr] : coordinate === 'y' ? [0, i, ddr] : [ddr, 0, i]))
                    },
                    {
                        p1: new Point3d(...(coordinate === 'x' ? [-i, 0, -ddr] : coordinate === 'y' ? [0, -i, -ddr] : [-ddr, 0, -i])),
                        p2: new Point3d(...(coordinate === 'x' ? [-i, 0, ddr] : coordinate === 'y' ? [0, -i, ddr] : [ddr, 0, -i]))
                    }
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

//  Drawing the observer's viewing frustum
//  Four rays are drawn from the observer position through the corners
//  of the observer's virtual screen.
//  First, the screen corners are described in the observer's local system
//  and transformed back to the global system using d_obs.
//  Then, those rays are projected onto your screen using d.
function draw_view_frustum(d, d_obs, observer) {
    // Screen corners in the observer's local coordinate system
    const p_screen = [
        new Point3d(-e_width / 2,  e_height / 2, e_dist),
        new Point3d( e_width / 2,  e_height / 2, e_dist),
        new Point3d( e_width / 2, -e_height / 2, e_dist),
        new Point3d(-e_width / 2, -e_height / 2, e_dist)
    ];

    const observerScreenPoint = d.point_3d(observer);
    if (!observerScreenPoint) {
        return;
    }

    const frustumColor = "rgb(255, 0, 0)";

    for (const corner of p_screen) {
        // Transform screen corner from observer space back to world space
        const cornerWorld = d_obs.inverse_transform(corner);

        // Direction from observer to that corner
        const ray = new Vector(observer, cornerWorld).multiplyByNumber(R_obs);

        // Extend the ray to get a visible frustum edge
        const farPoint = Point3d.translate(observer, ray);

        const farScreenPoint = d.point_3d(farPoint);
        if (!farScreenPoint) {
            return;
        }

        drawLine(
            observerScreenPoint.x, observerScreenPoint.y,
            farScreenPoint.x, farScreenPoint.y,
            frustumColor
        );
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

function draw_cube(d) {
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

function size_change(d) {
    _("size_info").innerHTML = d;
    w = d + "px";
    _(canvas_id).style.width = w;
    _(canvas_id).style.height = 3 * d / 4 + "px";
}

function set_zoom() {
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

function draw_scene() {
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
    draw_view_frustum(d, d_obs, observer);    // Draw the observer's viewing frustum edges
    draw_observer_screen(d, d_obs);           // Draw the observer's screen
    draw_cube(d);                             // Draw the cube
}
