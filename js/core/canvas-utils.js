function setPixel(x, y, r, g, b) {
    if (typeof g !== "undefined" && typeof b !== "undefined")
        ctx.fillStyle = "rgb(" + r + "," + g + "," + b + ")";
    else
        ctx.fillStyle = r;

    ctx.fillRect(x, y, 1, 1);
}

function drawLine(x1, y1, x2, y2, r, g, b) {
    if (typeof g !== "undefined" && typeof b !== "undefined")
        ctx.strokeStyle = "rgb(" + r + "," + g + "," + b + ")";
    else
        ctx.strokeStyle = r;

    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.lineWidth = 1;
    ctx.stroke();
}
