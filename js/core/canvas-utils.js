function setPixel(x, y, r, g, b) {
    if (typeof g !== "undefined" && typeof b !== "undefined")
        ctx.fillStyle = "rgb(" + r + "," + g + "," + b + ")";
    else
        ctx.fillStyle = r;

    ctx.fillRect(x, y, 1, 1);
}

function drawLineCanvas(x1, y1, x2, y2, r, g, b) {
    if (typeof g !== "undefined" && typeof b !== "undefined")
        ctx.strokeStyle = "rgb(" + r + "," + g + "," + b + ")";
    else
        ctx.strokeStyle = r;

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

function setPixelFast(x, y, r, g, b) {
    if (x < 0 || y < 0 || x >= width || y >= height)
        return;

    const i = (y * width + x) * 4;

    buffer[i] = r;
    buffer[i + 1] = g;
    buffer[i + 2] = b;
    buffer[i + 3] = 255;
}

function parseColor(c) {
    if (typeof c === "string") {
        const m = c.match(/\d+/g);
        return [parseInt(m[0]), parseInt(m[1]), parseInt(m[2])];
    }
    return c;
}

function drawLine(x0, y0, x1, y1, r, g, b) {
    if (typeof g === "undefined") {
        [r, g, b] = parseColor(r);
    }

    x0 = Math.round(x0);
    y0 = Math.round(y0);
    x1 = Math.round(x1);
    y1 = Math.round(y1);

    let dx = Math.abs(x1 - x0);
    let dy = Math.abs(y1 - y0);

    let sx = x0 < x1 ? 1 : -1;
    let sy = y0 < y1 ? 1 : -1;

    let err = dx - dy;

    while (true) {

        setPixelFast(x0, y0, r, g, b);

        if (x0 === x1 && y0 === y1)
            break;

        let e2 = 2 * err;

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
