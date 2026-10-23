// Render index.html (overview page) inside a minimal publish-like shell and save screenshots (light/dark x 1280/400 px).
// Reports console errors, horizontal overflow, and overlapping or clipped SVG text labels.
// Adapted from hyperloop/infographic/shoot.js (single page, no version tabs).
// Usage: PLAYWRIGHT_BROWSERS_PATH=/opt/pw-browsers NODE_PATH=/opt/node22/lib/node_modules node shoot.js [outDir]
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const here = __dirname;
const out = process.argv[2] || path.join(here, 'preview');
const page = fs.readFileSync(path.join(here, 'index.html'), 'utf8');
const shell = `<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top);padding-bottom:env(safe-area-inset-bottom)}
body{margin:0;font:14px system-ui}img{max-width:100%}[hidden]{display:none!important}</style>
</head><body>${page}</body></html>`;
const tmp = path.join(out, '_shell.html');
fs.mkdirSync(out, { recursive: true });
fs.writeFileSync(tmp, shell);

function audit() {
  const res = [];
  const scrolls = el => { for (let a = el.parentElement; a && a !== document.body; a = a.parentElement) { const o = getComputedStyle(a).overflowX; if (o === 'auto' || o === 'scroll' || o === 'hidden') return true; } return false; };
  document.querySelectorAll('svg').forEach((svg, si) => {
    const ts = [...svg.querySelectorAll('text')].map(t => ({ t: t.textContent, r: t.getBoundingClientRect() })).filter(o => o.r.width > 0);
    for (let i = 0; i < ts.length; i++) for (let j = i + 1; j < ts.length; j++) {
      const a = ts[i].r, b = ts[j].r;
      const ox = Math.min(a.right, b.right) - Math.max(a.left, b.left), oy = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
      if (ox > 1.5 && oy > 0.35 * Math.min(a.height, b.height)) res.push(`svg${si}: "${ts[i].t}" x "${ts[j].t}" (${ox.toFixed(0)}x${oy.toFixed(0)})`);
    }
    const sr = svg.getBoundingClientRect();
    ts.forEach(o => { if (o.r.left < sr.left - 1 || o.r.right > sr.right + 1 || o.r.top < sr.top - 1 || o.r.bottom > sr.bottom + 1) res.push(`svg${si}: "${o.t}" clipped`); });
    // smallest rendered font size
    const min = Math.min(...[...svg.querySelectorAll('text')].map(t => parseFloat(t.getAttribute('font-size')) * sr.width / svg.viewBox.baseVal.width));
    res.push(`svg${si}: min rendered font ${min.toFixed(1)} px`);
  });
  document.querySelectorAll('.wrap *').forEach(el => {
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.right > document.documentElement.clientWidth + 1 && !scrolls(el)) res.push('overflow: ' + el.tagName + '.' + el.className + ' ' + Math.round(r.right));
  });
  return res.slice(0, 40);
}

(async () => {
  const browser = await chromium.launch();
  let fails = 0;
  for (const scheme of ['light', 'dark']) {
    for (const width of [1280, 400]) {
      const ctx = await browser.newContext({ viewport: { width, height: 900 }, colorScheme: scheme, reducedMotion: 'reduce' });
      const p = await ctx.newPage();
      const errs = [];
      p.on('pageerror', e => errs.push(e.message));
      p.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
      await p.goto('file://' + tmp, { waitUntil: 'load' });
      await p.waitForTimeout(1200);
      const sw = await p.evaluate(() => [document.documentElement.scrollWidth, window.innerWidth]);
      const issues = await p.evaluate(audit);
      const file = path.join(out, `${scheme}-${width}.png`);
      await p.screenshot({ path: file, fullPage: true });
      console.log(path.basename(file), 'scrollWidth/innerWidth', sw.join('/'), errs.length ? 'ERRORS: ' + errs.join(' | ') : '');
      if (sw[0] > sw[1]) { fails++; console.log('  FAIL horizontal scroll'); }
      if (errs.length) fails++;
      issues.forEach(i => console.log('   ', i));
      await ctx.close();
    }
  }
  await browser.close();
  fs.unlinkSync(tmp);
  console.log(fails ? `${fails} failure(s)` : 'all checks passed');
})();
