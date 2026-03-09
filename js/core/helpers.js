function relativeDifference(a, b) {
    let absA = Math.abs(a);
    let absB = Math.abs(b);
    let maxAbs = Math.max(absA, absB);
    return maxAbs === 0 ? 0 : Math.abs(a - b) / maxAbs;
}

function sphericalToCartesian(r, phi, theta) {
    const k = Math.PI / 180.0;

    const x = r * Math.sin(k * theta) * Math.cos(k * phi);
    const y = r * Math.sin(k * theta) * Math.sin(k * phi);
    const z = r * Math.cos(k * theta);

    return new Point3d(x, y, z);
}

function _(id) {
    return document.getElementById(id);
}
