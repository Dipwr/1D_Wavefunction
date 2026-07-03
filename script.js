function simulator(maxX, numPoints) {
    this.maxX = maxX;
    this.numPoints = numPoints;
    this.dx = this.maxX / this.numPoints;
    this.psiRe = new Float64Array(numPoints);
    this.psiIm = new Float64Array(numPoints);
    this.hbar = 1;
    this.mass = 1;

    this.Ve = function(x) {
        return 0; 
    }

    this.initWavefuntion = function(x0, k0, sigma) {
        let sumProb = 0;
        for (let i = 0; i < this.numPoints; i++) {
            let x = i * this.dx;
            let envelope = Math.exp(-0.5 * Math.pow((x - x0) / sigma, 2));
            let phase = k0 * x;
            this.psiRe[i] = envelope * Math.cos(phase);
            this.psiIm[i] = envelope * Math.sin(phase);
            sumProb += (this.psiRe[i] * this.psiRe[i] + this.psiIm[i] * this.psiIm[i]) * this.dx;
        }

        let normFactor = Math.sqrt(sumProb);
        for (let i = 0; i < this.numPoints; i++) {
            this.psiRe[i] /= normFactor;
            this.psiIm[i] /= normFactor;
        }
    }

    this.step = function(dt) {
        const c = -(this.hbar * this.hbar) / (2 * this.mass * this.dx * this.dx);
        for (let i = 1; i < this.numPoints - 1; i++) {
            let x = i * this.dx;
            let d2psiIm = this.psiIm[i + 1] - 2 * this.psiIm[i] + this.psiIm[i - 1];
            let dpsiRe = (c * d2psiIm + this.Ve(x) * this.psiIm[i]) / this.hbar;
            this.psiRe[i] += dpsiRe * dt;
        }
        for (let i = 1; i < this.numPoints - 1; i++) {
            let x = i * this.dx;
            let d2psiRe = this.psiRe[i + 1] - 2 * this.psiRe[i] + this.psiRe[i - 1];
            let dpsiIm = -(c * d2psiRe + this.Ve(x) * this.psiRe[i]) / this.hbar;
            this.psiIm[i] += dpsiIm * dt;
        }
    }
}

const renderer = new THREE.WebGLRenderer({ canvas: document.getElementById("simulation"), antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(-40, 25, 50);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.target.set(0, 0, 0);

const simGroup = new THREE.Group();
const labelsGroup = new THREE.Group();
scene.add(simGroup);
scene.add(labelsGroup);

let sim, waveGeo, probGeo;
let arrows = [];
let waveScale = 15;
let probScale = 60;
let dt = 0.001;
let isPlaying = false;
let stepsPerFrame = 1;
let xOffset = 0;

function createTextSprite(text, fontSize, color) {
    const cvs = document.createElement('canvas');
    const ctx = cvs.getContext('2d');
    ctx.font = `bold ${fontSize}px monospace`;
    const metrics = ctx.measureText(text);
    cvs.width = Math.max(64, metrics.width + 20);
    cvs.height = fontSize + 20;
    
    ctx.font = `bold ${fontSize}px monospace`;
    ctx.fillStyle = color;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, cvs.width / 2, cvs.height / 2);
    
    const texture = new THREE.CanvasTexture(cvs);
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture }));
    sprite.scale.set(cvs.width / 15, cvs.height / 15, 1);
    return sprite;
}

function clearGroup(group) {
    while (group.children.length > 0) {
        const child = group.children[0];
        group.remove(child);
        if (child.geometry) child.geometry.dispose();
        if (child.material) {
            if (child.material.map) child.material.map.dispose();
            child.material.dispose();
        }
        if (child.children) clearGroup(child); 
    }
}

function rebuildAxes() {
    clearGroup(labelsGroup);

    for (let i = 0; i <= sim.maxX; i += 20) {
        const sprite = createTextSprite(i.toString(), 24, "#cccccc");
        sprite.position.set(i - xOffset, -3, 0);
        labelsGroup.add(sprite);
    }

    const visualSteps = [-30, -15, 15, 30];
    for (let vy of visualSteps) {
        const spriteRe = createTextSprite((vy / waveScale).toFixed(2), 20, "#55ff55");
        spriteRe.position.set(-xOffset - 4, vy, 0);
        labelsGroup.add(spriteRe);

        const spriteIm = createTextSprite((vy / waveScale).toFixed(2), 20, "#5555ff");
        spriteIm.position.set(-xOffset - 4, 0, vy);
        labelsGroup.add(spriteIm);
    }

    const probSteps = [15, 30, 45, 60];
    for (let vy of probSteps) {
        const spritePr = createTextSprite((vy / probScale).toFixed(3), 20, "#ffaa00");
        spritePr.position.set(-xOffset + 4, vy, 0);
        labelsGroup.add(spritePr);
    }
}

function buildSimulation() {
    clearGroup(simGroup);
	arrows = [];

	const maxX = parseFloat(document.getElementById("maxXSlider").value);
	const numPoints = parseInt(document.getElementById("ptsValue").innerText);
	const x0 = parseFloat(document.getElementById("posSlider").value);
	const k0 = parseFloat(document.getElementById("momSlider").value);
	const sigma = parseFloat(document.getElementById("sigmaSlider").value); 

	sim = new simulator(maxX, numPoints);
	sim.initWavefuntion(x0, k0, sigma); 
	xOffset = sim.maxX / 2;

    const axisMat = new THREE.LineBasicMaterial({ color: 0x444444 });
    simGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-xOffset, 0, 0), new THREE.Vector3(xOffset, 0, 0)
    ]), axisMat));
    simGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-xOffset, -60, 0), new THREE.Vector3(-xOffset, 60, 0)
    ]), axisMat));
    simGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-xOffset, 0, -60), new THREE.Vector3(-xOffset, 0, 60)
    ]), axisMat));

    const titleX = createTextSprite("Position (X)", 32, "#ff5555");
    titleX.position.set(xOffset + 5, 0, 0);
    simGroup.add(titleX);
    const titleY1 = createTextSprite("Amplitude (Y)", 32, "#55ff55");
    titleY1.position.set(-xOffset - 8, 40, 0);
    simGroup.add(titleY1);
    const titleY2 = createTextSprite("Probability (Y)", 32, "#ffaa00");
    titleY2.position.set(-xOffset + 8, 40, 0);
    simGroup.add(titleY2);
    const titleZ = createTextSprite("Imaginary (Z)", 32, "#5555ff");
    titleZ.position.set(-xOffset, 0, 40);
    simGroup.add(titleZ);

    waveGeo = new THREE.BufferGeometry();
    waveGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(sim.numPoints * 3), 3));
    simGroup.add(new THREE.Line(waveGeo, new THREE.LineBasicMaterial({ color: 0x00ffff, linewidth: 2 })));

    probGeo = new THREE.BufferGeometry();
    probGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(sim.numPoints * 3), 3));
    simGroup.add(new THREE.Line(probGeo, new THREE.LineBasicMaterial({ color: 0xffaa00, transparent: true, opacity: 0.8 })));

    const arrowSpacing = Math.max(1, Math.floor(sim.numPoints / 32)); 
    for (let i = 0; i < sim.numPoints; i += arrowSpacing) {
        const arrow = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3((i * sim.dx) - xOffset, 0, 0), 1, 0xaaaaaa, 0.5, 0.5);
        simGroup.add(arrow);
        arrows.push({ index: i, helper: arrow });
    }

    rebuildAxes();
    updateRenderData();
}

// Standard linear binder
function bindSlider(id, valueId, callback) {
    document.getElementById(id).addEventListener("input", (e) => {
        document.getElementById(valueId).innerText = e.target.value;
        if(callback) callback(parseFloat(e.target.value));
    });
}

// Logarithmic binder
function bindLogSlider(id, valueId, callback, isInt = false, fractionDigits = 2) {
    const slider = document.getElementById(id);
    const minVal = parseFloat(slider.getAttribute("data-min"));
    const maxVal = parseFloat(slider.getAttribute("data-max"));
    
    slider.addEventListener("input", (e) => {
        // Normalize the 0-1000 slider to a 0.0 - 1.0 percentage
        const t = parseFloat(e.target.value) / 1000;
        
        // Exponential mapping formula: val = min * (max/min)^t
        let actualVal = minVal * Math.pow(maxVal / minVal, t);
        
        if (isInt) {
            actualVal = Math.round(actualVal);
            document.getElementById(valueId).innerText = actualVal;
        } else {
            document.getElementById(valueId).innerText = actualVal.toFixed(fractionDigits);
        }
        
        if(callback) callback(actualVal);
    });
}

document.getElementById("playPauseBtn").addEventListener("click", (e) => {
    isPlaying = !isPlaying;
    e.target.innerText = isPlaying ? "Pause" : "Play";
    document.getElementById("stepBtn").disabled = isPlaying;
});

document.getElementById("stepBtn").addEventListener("click", () => {
    if (!isPlaying) { sim.step(dt); updateRenderData(); }
});

document.getElementById("resetBtn").addEventListener("click", () => {
    buildSimulation();
});

// --- INITIALIZE LINEAR SLIDERS ---
bindSlider("posSlider", "posValue");
bindSlider("momSlider", "momValue");
bindSlider("sigmaSlider", "sigmaValue");
bindSlider("maxXSlider", "maxXValue");

// --- INITIALIZE LOGARITHMIC SLIDERS ---
// dt uses 4 decimal places
bindLogSlider("dtSlider", "dtValue", (val) => dt = val, false, 4); 

// Speed uses integers
bindLogSlider("speedSlider", "speedValue", (val) => stepsPerFrame = val, true); 

// Resolution uses integers
bindLogSlider("ptsSlider", "ptsValue", null, true); 

// Scales use integers and rebuild the axes dynamically
bindLogSlider("wsSlider", "wsValue", (val) => {
    waveScale = val;
    rebuildAxes();
    if (!isPlaying) updateRenderData();
}, true);

bindLogSlider("psSlider", "psValue", (val) => {
    probScale = val;
    rebuildAxes();
    if (!isPlaying) updateRenderData();
}, true);

window.addEventListener('resize', () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
});

function updateRenderData() {
    const wPos = waveGeo.attributes.position.array;
    const pPos = probGeo.attributes.position.array;

    for (let i = 0; i < sim.numPoints; i++) {
        const x3D = (i * sim.dx) - xOffset;
        const re = sim.psiRe[i];
        const im = sim.psiIm[i];
        const prob = (re * re) + (im * im);

        const idx = i * 3;
        wPos[idx] = x3D; wPos[idx + 1] = re * waveScale; wPos[idx + 2] = im * waveScale;
        pPos[idx] = x3D; pPos[idx + 1] = prob * probScale; pPos[idx + 2] = 0; 
    }
    
    waveGeo.attributes.position.needsUpdate = true;
    probGeo.attributes.position.needsUpdate = true;

    for (const a of arrows) {
        const re = sim.psiRe[a.index];
        const im = sim.psiIm[a.index];
        let mag = Math.sqrt(re*re + im*im) * waveScale;
        if (mag < 0.01) mag = 0.01; 
        a.helper.setDirection(new THREE.Vector3(0, re, im).normalize());
        a.helper.setLength(mag, mag * 0.2, mag * 0.2);
    }
}

function frame() {
    if (isPlaying) {
        for (let i = 0; i < stepsPerFrame; i++) sim.step(dt);
        updateRenderData();
    }
    controls.update();
    renderer.render(scene, camera);
    requestAnimationFrame(frame);
}

buildSimulation();
frame();