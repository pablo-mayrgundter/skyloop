// Render index.html (BoomShield infographic) inside a minimal publish-like shell and save screenshots
// for every concept version tab (light/dark x 1280/400 px). Also reports console errors, horizontal
// overflow, overlapping or clipped SVG text labels, and checks the version tabs (click, keyboard, hash).
// Same script on the BoomShield, Skyloop and Hyperloop pages.
// Usage: PLAYWRIGHT_BROWSERS_PATH=/opt/pw-browsers NODE_PATH=/opt/node22/lib/node_modules node shoot.js [outDir]
// Output: preview/<scheme>-<width>.png for the default (latest) version, preview/<scheme>-<width>-<version>.png for earlier ones.
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const here = __dirname;
const out = process.argv[2] || path.join(here, 'preview');
const page = fs.readFileSync(path.join(here, 'index.html'), 'utf8');
const shell = `<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top);padding-bottom:env(safe-area-inset-bottom)}
body{margin:0;font:14px system-ui}img{max-width:100%}</style>
</head><body>${page}</body></html>`;
const tmp = path.join(out, '_shell.html');
fs.mkdirSync(out, { recursive: true });
fs.writeFileSync(tmp, shell);
const url = 'file://' + tmp;

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
    ts.forEach(o => { if (o.r.left < sr.left - 1 || o.r.right > sr.right + 1) res.push(`svg${si}: "${o.t}" clipped`); });
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
  const bad = m => { fails++; console.log('  FAIL', m); };

  // versions offered by the tab bar (first = latest = default)
  const probe = await (await browser.newContext()).newPage();
  await probe.goto(url, { waitUntil: 'load' });
  const versions = await probe.$$eval('#vtabs [role="tab"]', ts => ts.map(t => t.dataset.v));
  console.log('versions:', versions.join(', '));

  for (const scheme of ['light', 'dark']) {
    for (const width of [1280, 400]) {
      for (const [vi, ver] of versions.entries()) {
        const ctx = await browser.newContext({ viewport: { width, height: 900 }, colorScheme: scheme, reducedMotion: 'reduce' });
        const p = await ctx.newPage();
        const errs = [];
        p.on('pageerror', e => errs.push(e.message));
        p.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
        await p.goto(url + (vi ? '#' + ver : ''), { waitUntil: 'load' });
        await p.waitForTimeout(1000);
        const st = await p.evaluate(() => ({
          sw: document.documentElement.scrollWidth, iw: window.innerWidth,
          sel: document.querySelector('#vtabs [aria-selected="true"]')?.dataset.v,
          shown: [...document.querySelectorAll('[role="tabpanel"]')].filter(x => !x.hidden).map(x => x.id)
        }));
        const issues = await p.evaluate(audit);
        const file = path.join(out, `${scheme}-${width}${vi ? '-' + ver : ''}.png`);
        await p.screenshot({ path: file, fullPage: true });
        console.log(path.basename(file), 'scrollWidth/innerWidth', st.sw + '/' + st.iw, 'tab', st.sel, 'panel', st.shown.join(','), errs.length ? 'ERRORS: ' + errs.join(' | ') : '');
        if (st.sw > st.iw) bad('horizontal scroll');
        if (st.sel !== ver || st.shown.length !== 1) bad('wrong tab/panel for ' + ver);
        if (errs.length) fails++;
        issues.forEach(i => console.log('   ', i));
        await ctx.close();
      }
    }
  }

  // interaction checks: click, keyboard, hashchange, in-page version links, redraw after a hidden render
  for (const width of [1280, 400]) {
    const ctx = await browser.newContext({ viewport: { width, height: 900 }, reducedMotion: 'reduce' });
    const p = await ctx.newPage();
    const errs = []; p.on('pageerror', e => errs.push(e.message));
    const state = () => p.evaluate(() => ({
      sel: document.querySelector('#vtabs [aria-selected="true"]')?.dataset.v,
      shown: [...document.querySelectorAll('[role="tabpanel"]')].filter(x => !x.hidden).map(x => x.id).join(','),
      hash: location.hash, eyebrow: document.querySelector('header .eyebrow').textContent.trim(),
      focus: document.activeElement?.dataset?.v || null, heroW: document.querySelector('#hero svg')?.viewBox.baseVal.width || null
    }));
    await p.goto(url, { waitUntil: 'load' });
    let s = await state();
    console.log(`[${width}] default:`, JSON.stringify(s));
    if (s.sel !== versions[0]) bad('default is not the latest version');
    if (versions.length > 1) {
      await p.click(`#vtabs [data-v="${versions[1]}"]`); s = await state();
      console.log(`[${width}] click ${versions[1]}:`, JSON.stringify(s));
      if (s.sel !== versions[1] || s.hash !== '#' + versions[1] || !s.eyebrow.includes(versions[1])) bad('click did not switch tab/hash/eyebrow');
      await p.keyboard.press('ArrowLeft'); s = await state();
      console.log(`[${width}] ArrowLeft:`, JSON.stringify(s));
      if (s.sel !== versions[0] || s.focus !== versions[0]) bad('keyboard did not move selection/focus');
      await p.keyboard.press('End'); s = await state();
      if (s.sel !== versions[versions.length - 1]) bad('End key');
      await p.keyboard.press('Home'); s = await state();
      if (s.sel !== versions[0]) bad('Home key');
      await p.evaluate(v => { location.hash = '#' + v; }, versions[1]); await p.waitForTimeout(100); s = await state();
      console.log(`[${width}] hashchange ${versions[1]}:`, JSON.stringify(s));
      if (s.sel !== versions[1]) bad('hashchange did not switch');
      // in-page "back to latest" link, then "see earlier" link
      await p.click(`[role="tabpanel"]:not([hidden]) a.vgo`); s = await state();
      if (s.sel !== versions[0]) bad('vgo link (back) did not switch');
      await p.click(`[role="tabpanel"]:not([hidden]) a.vgo`); s = await state();
      console.log(`[${width}] vgo links:`, JSON.stringify(s));
      if (s.sel !== versions[1]) bad('vgo link (earlier) did not switch');
      // deep link to the earlier version, then show the latest: width-dependent figures must be redrawn
      const q = await ctx.newPage(); q.on('pageerror', e => errs.push(e.message));
      await q.goto(url + '#' + versions[1], { waitUntil: 'load' });
      await q.click(`#vtabs [data-v="${versions[0]}"]`);
      const hw = await q.evaluate(() => document.querySelector('#hero svg')?.viewBox.baseVal.width || null);
      console.log(`[${width}] deep link ${versions[1]} -> ${versions[0]}: hero viewBox width`, hw, 'vs', s.heroW);
      await q.close();
    } else {
      await p.keyboard.press('Tab');
      await p.evaluate(v => { location.hash = '#' + v; }, versions[0]); s = await state();
      if (s.sel !== versions[0]) bad('single-version hash');
    }
    if (errs.length) bad('page errors: ' + errs.join(' | '));
    await ctx.close();
  }
  await browser.close();
  fs.unlinkSync(tmp);
  console.log(fails ? `${fails} failure(s)` : 'all checks passed');
})();
