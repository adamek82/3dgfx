const COLOR_BLACK      = [0, 0, 0];
const COLOR_RED        = [255, 0, 0];
const COLOR_GREEN      = [0, 255, 0];
const COLOR_BLUE       = [0, 0, 255];
const COLOR_YELLOW     = [255, 255, 0];
const COLOR_LIGHT_GRAY = [192, 192, 192];

function colorToCss(color) {
    return `rgb(${color[0]}, ${color[1]}, ${color[2]})`;
}

function setPixel(x, y, color) {
    ctx.fillStyle = colorToCss(color);
    ctx.fillRect(x, y, 1, 1);
}

function drawLineCanvas(x1, y1, x2, y2, color) {
    ctx.strokeStyle = colorToCss(color);

    ctx.beginPath();
    ctx.moveTo(x1 + 0.5, y1 + 0.5);
    ctx.lineTo(x2 + 0.5, y2 + 0.5);
    ctx.stroke();
}

function initPixelBuffer(w, h, imageData) {
    width = w;
    height = h;
    buffer = imageData.data;
}

function setPixelFast(x, y, color) {
    if (x < 0 || y < 0 || x >= width || y >= height)
        return;

    const i = (y * width + x) * 4;

    buffer[i]     = color[0];
    buffer[i + 1] = color[1];
    buffer[i + 2] = color[2];
    buffer[i + 3] = 255;
}

function parseColor(c) {
    if (typeof c === "string") {
        const m = c.match(/\d+/g);
        return [parseInt(m[0]), parseInt(m[1]), parseInt(m[2])];
    }
    return c;
}

function drawLine(x0, y0, x1, y1, color) {
    x0 |= 0;
    y0 |= 0;
    x1 |= 0;
    y1 |= 0;

    let dx = x1 - x0;
    let dy = y1 - y0;

    const sx = dx >= 0 ? 1 : -1;
    const sy = dy >= 0 ? 1 : -1;

    dx = dx * sx;
    dy = dy * sy;

    let err = dx - dy;

    while (true) {
        setPixelFast(x0, y0, color);

        if (x0 === x1 && y0 === y1)
            break;

        const e2 = err << 1;   // faster than 2 * err

        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }

        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
    }
}
