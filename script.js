function simulator(maxX, numPoints) {
	this.maxX = maxX;
	this.numPoints = numPoints;
	this.dx = this.maxX / this.numPoints;
	this.psiRe = new Float64Array(numPoints);
	this.psiIm = new Float64Array(numPoints);
	
	this.cPrimeRe = new Float64Array(numPoints);
	this.cPrimeIm = new Float64Array(numPoints);
	this.dPrimeRe = new Float64Array(numPoints);
	this.dPrimeIm = new Float64Array(numPoints);
	
	this.hbar = 1;
	this.mass = 1;
	this.regions = []; // Master list of potential regions

	this.Ve = function(x) {
		for (let i = 0; i < this.regions.length; i++) {
			let r = this.regions[i];
			if (x >= r.startX && x <= r.endX) {
				if (r.type === 'constant') return r.params.height;
				if (r.type === 'linear') return r.params.slope * (x - r.startX) + r.params.base;
				if (r.type === 'harmonic') return 0.5 * r.params.k * Math.pow(x - r.params.center, 2);
				if (r.type === 'custom') return r._compiledFn ? r._compiledFn(x) : 0;
			}
		}
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
		const r = (this.hbar * dt) / (4 * this.mass * this.dx * this.dx);

		for (let i = 1; i < this.numPoints - 1; i++) {
			let vj = (dt / (2 * this.hbar)) * this.Ve(i * this.dx);

			let B_Re = 1.0;
			let B_Im = 2 * r + vj;

			let old_j_Re = this.psiRe[i];
			let old_j_Im = this.psiIm[i];
			let old_jm1_Re = this.psiRe[i-1];
			let old_jm1_Im = this.psiIm[i-1];
			let old_jp1_Re = this.psiRe[i+1];
			let old_jp1_Im = this.psiIm[i+1];

			let term2_Re = old_j_Re - (-(2 * r + vj)) * old_j_Im;
			let term2_Im = old_j_Im + (-(2 * r + vj)) * old_j_Re;

			let term1_Re = -r * old_jm1_Im;
			let term1_Im =  r * old_jm1_Re;
			let term3_Re = -r * old_jp1_Im;
			let term3_Im =  r * old_jp1_Re;

			let D_Re = term1_Re + term2_Re + term3_Re;
			let D_Im = term1_Im + term2_Im + term3_Im;

			if (i === 1) {
				let denom = B_Re * B_Re + B_Im * B_Im;
				this.cPrimeRe[i] = (-r * B_Im) / denom;
				this.cPrimeIm[i] = (-r * B_Re) / denom;
				this.dPrimeRe[i] = (D_Re * B_Re + D_Im * B_Im) / denom;
				this.dPrimeIm[i] = (D_Im * B_Re - D_Re * B_Im) / denom;
			} else {
				let sub_Re =  r * this.cPrimeIm[i-1];
				let sub_Im = -r * this.cPrimeRe[i-1];
				
				let den_Re = B_Re - sub_Re;
				let den_Im = B_Im - sub_Im;
				let denomMag = den_Re * den_Re + den_Im * den_Im;

				this.cPrimeRe[i] = (-r * den_Im) / denomMag;
				this.cPrimeIm[i] = (-r * den_Re) / denomMag;

				let num_Re = D_Re - (r * this.dPrimeIm[i-1]);
				let num_Im = D_Im - (-r * this.dPrimeRe[i-1]);

				this.dPrimeRe[i] = (num_Re * den_Re + num_Im * den_Im) / denomMag;
				this.dPrimeIm[i] = (num_Im * den_Re - num_Re * den_Im) / denomMag;
			}
		}

		for (let i = this.numPoints - 2; i >= 1; i--) {
			let next_Re = this.psiRe[i+1];
			let next_Im = this.psiIm[i+1];

			let mult_Re = this.cPrimeRe[i] * next_Re - this.cPrimeIm[i] * next_Im;
			let mult_Im = this.cPrimeRe[i] * next_Im + this.cPrimeIm[i] * next_Re;

			this.psiRe[i] = this.dPrimeRe[i] - mult_Re;
			this.psiIm[i] = this.dPrimeIm[i] - mult_Im;
		}
	}

	this.calculateEnergies = function() {
		let expV = 0;
		let expK = 0;
		const C = (this.hbar * this.hbar) / (2 * this.mass * this.dx * this.dx);

		for (let i = 1; i < this.numPoints - 1; i++) {
			let re = this.psiRe[i];
			let im = this.psiIm[i];
			let prob = re * re + im * im;

			// Expected Potential <V> = Integral( P(x) * V(x) dx )
			expV += prob * this.Ve(i * this.dx) * this.dx;

			// Expected Kinetic <K> = Integral( Psi* (-hbar^2/2m d^2/dx^2) Psi dx )
			let d2re = this.psiRe[i+1] - 2 * re + this.psiRe[i-1];
			let d2im = this.psiIm[i+1] - 2 * im + this.psiIm[i-1];
			
			let kRe = -C * d2re;
			let kIm = -C * d2im;

			// Multiply by complex conjugate (Re - i*Im)
			expK += (re * kRe + im * kIm) * this.dx;
		}
		
		return { k: expK, v: expV, t: expK + expV };
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

let sim, waveGeo, probGeo, potGeo;
let arrows = [];
let waveScale = 15;
let probScale = 60;
let potScale = 2.0;
let dt = 0.01;
let isPlaying = false;
let isRealtime = true;
let lastFrameTime = 0;
let stepsPerFrame = 1;
let subSteps = 25;
let xOffset = 0;

// Default initial barrier
let potentialRegions = [
	{ startX: 45, endX: 55, type: 'constant', params: { height: 10 } }
];

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

function updatePotentialVisuals() {
	if (!potGeo) return;
	const pos = potGeo.attributes.position.array;
	for (let i = 0; i < sim.numPoints; i++) {
		const x3D = (i * sim.dx) - xOffset;
		const v = sim.Ve(i * sim.dx);
		const idx = i * 3;
		pos[idx] = x3D;
		pos[idx + 1] = v * potScale; 
		pos[idx + 2] = -0.1; // Push slightly behind the wave
	}
	potGeo.attributes.position.needsUpdate = true;
	updateEnergyReadout();
}

function buildSimulation() {
	clearGroup(simGroup);
	arrows = [];

	const maxX = parseFloat(document.getElementById("maxXInput").value);
	const numPoints = parseInt(document.getElementById("ptsInput").value);
	const x0 = parseFloat(document.getElementById("posInput").value);
	const k0 = parseFloat(document.getElementById("momInput").value);
	const sigma = parseFloat(document.getElementById("sigmaInput").value); 

	sim = new simulator(maxX, numPoints);
	sim.regions = potentialRegions; // Sync UI array to engine
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
	
	// Add Potential Energy Title
	const titleV = createTextSprite("Potential V(x)", 32, "#ffffff");
	titleV.position.set(xOffset + 5, 20, 0);
	simGroup.add(titleV);

	waveGeo = new THREE.BufferGeometry();
	waveGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(sim.numPoints * 3), 3));
	simGroup.add(new THREE.Line(waveGeo, new THREE.LineBasicMaterial({ color: 0x00ffff, linewidth: 2 })));

	probGeo = new THREE.BufferGeometry();
	probGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(sim.numPoints * 3), 3));
	simGroup.add(new THREE.Line(probGeo, new THREE.LineBasicMaterial({ color: 0xffaa00, transparent: true, opacity: 0.8 })));

	// Instantiate the 3D potential line
	potGeo = new THREE.BufferGeometry();
	potGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(sim.numPoints * 3), 3));
	simGroup.add(new THREE.Line(potGeo, new THREE.LineBasicMaterial({ color: 0xffffff, linewidth: 2, transparent: true, opacity: 0.6 })));

	const arrowSpacing = Math.max(1, Math.floor(sim.numPoints / 32)); 
	for (let i = 0; i < sim.numPoints; i += arrowSpacing) {
		const arrow = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3((i * sim.dx) - xOffset, 0, 0), 1, 0xaaaaaa, 0.5, 0.5);
		simGroup.add(arrow);
		arrows.push({ index: i, helper: arrow });
	}

	rebuildAxes();
	updatePotentialVisuals();
	updateRenderData();
	updateEnergyReadout();
}

// --- POTENTIAL UI MANAGER ---
function renderRegionsUI() {
	const container = document.getElementById("regions-container");
	container.innerHTML = "";
	
	potentialRegions.forEach((reg, index) => {
		const card = document.createElement("div");
		card.className = "region-card";

		let paramsHTML = "";
		if (reg.type === 'constant') {
			paramsHTML = `<div><label>Height (V)</label><br><input type="number" step="0.5" value="${reg.params.height}" onchange="updateRegionParam(${index}, 'height', this.value)"></div>`;
		} else if (reg.type === 'linear') {
			paramsHTML = `
				<div><label>Base (V)</label><br><input type="number" step="0.5" value="${reg.params.base}" onchange="updateRegionParam(${index}, 'base', this.value)"></div>
				<div><label>Slope</label><br><input type="number" step="0.1" value="${reg.params.slope}" onchange="updateRegionParam(${index}, 'slope', this.value)"></div>
			`;
		} else if (reg.type === 'harmonic') {
			paramsHTML = `
				<div><label>Center X</label><br><input type="number" step="1" value="${reg.params.center}" onchange="updateRegionParam(${index}, 'center', this.value)"></div>
				<div><label>Spring (k)</label><br><input type="number" step="0.01" value="${reg.params.k}" onchange="updateRegionParam(${index}, 'k', this.value)"></div>
			`;
		} else if (reg.type === 'custom') {
			paramsHTML = `
				<div style="flex: 2;">
					<label>Formula: V(x) = </label><br>
					<input type="text" value="${reg.params.formula}" onchange="updateRegionParam(${index}, 'formula', this.value)" placeholder="e.g. 5 * Math.sin(x)">
					<div style="font-size: 10px; color: #888; margin-top: 4px;">Variables: 'x'. Available: Math.sin(), Math.cos(), Math.exp()</div>
				</div>
			`;
		}

		card.innerHTML = `
			<div class="region-header">
				Region ${index + 1}
				<button class="del-btn" onclick="deleteRegion(${index})">X</button>
			</div>
			<div class="region-row">
				<div><label>Start X</label><br><input type="number" value="${reg.startX}" step="1" onchange="updateRegion(${index}, 'startX', this.value)"></div>
				<div><label>End X</label><br><input type="number" value="${reg.endX}" step="1" onchange="updateRegion(${index}, 'endX', this.value)"></div>
			</div>
			<div class="region-row">
				<div style="flex: 2;"><label>Type</label><br>
					<select onchange="changeRegionType(${index}, this.value)">
						<option value="constant" ${reg.type === 'constant' ? 'selected' : ''}>Constant</option>
						<option value="linear" ${reg.type === 'linear' ? 'selected' : ''}>Linear Slope</option>
						<option value="harmonic" ${reg.type === 'harmonic' ? 'selected' : ''}>Harmonic Well</option>
						<option value="custom" ${reg.type === 'custom' ? 'selected' : ''}>Custom Formula</option>
					</select>
				</div>
			</div>
			<div class="region-row">${paramsHTML}</div>
		`;
		container.appendChild(card);
	});
}

window.compileCustomRegion = function(idx) {
	const reg = potentialRegions[idx];
	if (reg.type === 'custom') {
		try {
			// Create a fast, native JS function from the string
			reg._compiledFn = new Function('x', 'return ' + reg.params.formula + ';');
			// Test it immediately to catch syntax errors
			reg._compiledFn(0); 
		} catch (e) {
			console.warn("Invalid formula in Region " + (idx+1) + ". Defaulting to 0.");
			reg._compiledFn = function() { return 0; };
		}
	}
};

window.updateRegion = function(idx, field, val) {
	potentialRegions[idx][field] = parseFloat(val);
	updatePotentialVisuals();
};

window.updateRegionParam = function(idx, field, val) {
	if (field === 'formula') {
		potentialRegions[idx].params[field] = val;
		compileCustomRegion(idx); // Recompile on text change
	} else {
		potentialRegions[idx].params[field] = parseFloat(val);
	}
	updatePotentialVisuals();
};

window.changeRegionType = function(idx, type) {
	potentialRegions[idx].type = type;
	if (type === 'constant') potentialRegions[idx].params = { height: 10 };
	else if (type === 'linear') potentialRegions[idx].params = { base: 0, slope: 1 };
	else if (type === 'harmonic') potentialRegions[idx].params = { center: (potentialRegions[idx].startX + potentialRegions[idx].endX)/2, k: 0.1 };
	else if (type === 'custom') {
		potentialRegions[idx].params = { formula: "5 * Math.sin(x)" };
		compileCustomRegion(idx);
	}
	renderRegionsUI();
	updatePotentialVisuals();
};

window.deleteRegion = function(idx) {
	potentialRegions.splice(idx, 1);
	renderRegionsUI();
	updatePotentialVisuals();
};

document.getElementById("addRegionBtn").addEventListener("click", () => {
	potentialRegions.push({ startX: 40, endX: 60, type: 'constant', params: { height: 10 } });
	renderRegionsUI();
	updatePotentialVisuals();
});

// --- UI BINDINGS ---

function bindLogSlider(id, inputId, callback, isInt = false, fractionDigits = 4) {
	const slider = document.getElementById(id);
	const input = document.getElementById(inputId);
	const minVal = parseFloat(slider.getAttribute("data-min"));
	const maxVal = parseFloat(slider.getAttribute("data-max"));

	const valToSlider = (v) => Math.round(1000 * Math.log(v / minVal) / Math.log(maxVal / minVal));
	const sliderToVal = (s) => minVal * Math.pow(maxVal / minVal, s / 1000);

	const update = (val, fromSlider) => {
		if (fromSlider) {
			input.value = isInt ? Math.round(val) : val.toFixed(fractionDigits);
		} else {
			// Visually clamp the slider handle to the edges if they type outside bounds
			let clampedForSlider = Math.max(minVal, Math.min(maxVal, val));
			slider.value = valToSlider(clampedForSlider);
			
			// Format the input box but KEEP the out-of-bounds value for the physics engine
			input.value = isInt ? Math.round(val) : val.toFixed(fractionDigits);
		}
		if (callback) callback(val);
	};

	slider.addEventListener("input", (e) => update(sliderToVal(e.target.value), true));
	input.addEventListener("change", (e) => update(parseFloat(e.target.value), false));
	
	update(parseFloat(input.value), false);
}

function bindLinearSlider(id, inputId, callback) {
	const slider = document.getElementById(id);
	const input = document.getElementById(inputId);
	
	// Read the native min/max from the HTML element
	const minVal = parseFloat(slider.min);
	const maxVal = parseFloat(slider.max);

	const update = (val, fromSlider) => {
		if (fromSlider) {
			input.value = val;
		} else {
			// Visually clamp the slider handle to the edges
			slider.value = Math.max(minVal, Math.min(maxVal, val));
			input.value = val;
		}
		if (callback) callback(val);
	};

	slider.addEventListener("input", (e) => update(parseFloat(e.target.value), true));
	input.addEventListener("change", (e) => update(parseFloat(e.target.value), false));
	
	update(parseFloat(input.value), false);
}

renderRegionsUI();
buildSimulation();

bindLogSlider("dtSlider", "dtInput", (val) => dt = val, false, 4); 
bindLogSlider("wsSlider", "wsInput", (val) => { waveScale = val; rebuildAxes(); if(!isPlaying) updateRenderData(); }, true);
bindLogSlider("psSlider", "psInput", (val) => { probScale = val; rebuildAxes(); if(!isPlaying) updateRenderData(); }, true);
bindLogSlider("potSlider", "potInput", (val) => { potScale = val; updatePotentialVisuals(); }, false, 1);
bindLogSlider("speedSlider", "speedInput", (val) => stepsPerFrame = val, true);
bindLogSlider("ptsSlider", "ptsInput", null, true);
bindLogSlider("subSlider", "subInput", (val) => subSteps = Math.round(val), true);

bindLinearSlider("maxXSlider", "maxXInput", null, false);
bindLinearSlider("posSlider", "posInput", updateEnergyReadout, false);
bindLinearSlider("momSlider", "momInput", updateEnergyReadout, false);
bindLinearSlider("sigmaSlider", "sigmaInput", updateEnergyReadout, false);

document.getElementById("playPauseBtn").addEventListener("click", (e) => {
	isPlaying = !isPlaying;
	e.target.innerText = isPlaying ? "Pause" : "Play";
	document.getElementById("stepBtn").disabled = isPlaying;
	if (isPlaying && isRealtime) lastFrameTime = performance.now();
});

document.getElementById("stepBtn").addEventListener("click", () => {
	if (!isPlaying) { sim.step(dt); updateRenderData(); }
});

document.getElementById("resetBtn").addEventListener("click", () => {
	buildSimulation();
});

document.getElementById("realtimeCheck").addEventListener("change", (e) => {
	isRealtime = e.target.checked;
	document.getElementById("dtSlider").disabled = isRealtime;
	document.getElementById("dtInput").disabled = isRealtime;
	if (isRealtime && isPlaying) lastFrameTime = performance.now();
});

window.addEventListener('resize', () => {
	camera.aspect = window.innerWidth / window.innerHeight;
	camera.updateProjectionMatrix();
	renderer.setSize(window.innerWidth, window.innerHeight);
});

function updateEnergyReadout() {
	if (!sim) return;
	const energies = sim.calculateEnergies();

	document.getElementById("kin-eng").innerText = energies.k.toFixed(3);
	document.getElementById("pot-eng").innerText = energies.v.toFixed(3);
	document.getElementById("tot-eng").innerText = energies.t.toFixed(3);
}

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
	updateEnergyReadout();
}

let frameCount = 0;
let lastFpsUpdate = 0;

const fpsVal = document.getElementById("fps-val");
const msVal = document.getElementById("ms-val");

function frame() {
	let now = performance.now();
	let t0 = performance.now();

	if (isPlaying) {
		const totalTimePerFrame = (isRealtime ? (now - lastFrameTime) / 1000 : dt) * stepsPerFrame;
		const localDt = totalTimePerFrame / subSteps;
		
		for (let i = 0; i < subSteps; i++) {
			sim.step(localDt);
		}
		updateRenderData();
	}
	
	lastFrameTime = now;
	controls.update();
	renderer.render(scene, camera);
	
	let t1 = performance.now();
	msVal.innerText = (t1 - t0).toFixed(1);

	frameCount++;
	if (now - lastFpsUpdate >= 1000) {
		fpsVal.innerText = frameCount;
		frameCount = 0;
		lastFpsUpdate = now;
	}

	requestAnimationFrame(frame);
}

frame();