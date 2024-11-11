const MAXX = 1200;
const MAXY = 800;

const SCREEN_WIDTH = 0.3;
const SCREEN_HEIGHT = 0.2;
const DIST_SCREEN = 0.5;

var canvas_id = "mainCanvas",
	canvas = document.getElementById(canvas_id),
	ctx = canvas.getContext("2d"),
	width,
	height,
	X_obs = 200,
	Y_obs = 300,
	Z_obs = -30;
	points = [];

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

function point_through_window(x, y, z, X, Y, Z, color)
{
	A = width / SCREEN_WIDTH;
	B = width / 2;
	C = -height / SCREEN_HEIGHT;
	D = height / 2;
	
	r = Math.sqrt(X*X + Y*Y + Z*Z);
	v = Math.sqrt(X*X + Y*Y);
	
	if (v == 0) {
		sf = 0
		cf = 1
	} else {
		sf = Y / v;
		cf = X / v;
	}
	
	st = v / r;
	ct = Z / r;
	
	x_eye = -sf * x + cf * y;
	y_eye = -ct * (cf * x + sf * y) + st * z;
	z_eye = -st * (cf * x + sf * y) - ct * z + r;
	
	xe = A * x_eye * DIST_SCREEN / z_eye + B;
	ye = C * y_eye * DIST_SCREEN / z_eye + D;
	
	setPixel(xe, ye, color);
}

function draw_cube()
{
	ctx.fillStyle = 'black';
	ctx.fillRect(0, 0, MAXX, MAXY);
	
	let points = [
		[-5, -5, -5, 'blue'],
		[ 5, -5, -5, 'yellow'],
		[ 5,  5, -5, 'red'],
		[-5,  5, -5, 'blue'],
		[-5, -5,  5, 'green'],
		[ 5, -5,  5, 'yellow'],
		[ 5,  5,  5, 'red'],
		[-5,  5,  5, 'green'],
	];
	
	let X = 75; Y = 25; Z = 45;
	
	for (let i = 0; i < points.length; i++)
		point_through_window(points[i][0], points[i][1], points[i][2], X, Y, Z, points[i][3]);
}

function glass_skyscraper()
{
	if (points.length > 0)
		return;

	let no_of_points = 10000;
	let start_y = -170;
	let height = 60;

	for (let k = 0; k < 6; ++k) {
		let z = 20 + Math.random() * height;
		let wall_x = 15
		let wall_y = 25 + Math.floor(Math.random() * 25);
		for (let i = 0; i < Math.floor(no_of_points / 2); ++i) {
			points.push(
				[
					wall_x,
					Math.floor(Math.random() * wall_y) + start_y,
					Math.floor(Math.random() * z) - 40,
					'white'
				]
			);
		}
		for (let i = 0; i < Math.floor(no_of_points / 2); ++i) {
			points.push(
				[
					Math.floor(Math.random() * wall_x),
					wall_y + start_y,
					Math.floor(Math.random() * z) - 40,
					'gray'
				]
			);
		}
		start_y += wall_y + 15;
	}

}

function size_change(d)
{
	_("size_info").innerHTML = d;
	w = d + "px";
	_(canvas_id).style.width = w;
	_(canvas_id).style.height = 2 * d / 3 + "px";

}

function set_zoom()
{
	width = canvas.width = canvas.clientWidth;
	height = canvas.height = canvas.clientHeight;
}

function get_settings() {
    width = canvas.width;
    height = canvas.height;

	X_obs = _("X_obs").value;
    _("X_obs_info").innerHTML = X_obs;
	Y_obs = _("Y_obs").value;
    _("Y_obs_info").innerHTML = Y_obs;
	Z_obs = _("Z_obs").value;
    _("Z_obs_info").innerHTML = Z_obs;

    draw_scene()
}

function draw_scene()
{
	glass_skyscraper();
	
	set_zoom();
	
	ctx.fillStyle = 'black';
	ctx.fillRect(0, 0, width, height);

	for (let i = 0; i < points.length; i++)
		point_through_window(points[i][0], points[i][1], points[i][2], X_obs, Y_obs, Z_obs, points[i][3]);
}
