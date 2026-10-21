// Render index.html inside a minimal publish-like shell and save screenshots.
// Usage: PLAYWRIGHT_BROWSERS_PATH=/opt/pw-browsers node shoot.js [outDir]
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

(async () => {
  const browser = await chromium.launch();
  for (const scheme of ['light', 'dark']) {
    for (const width of [1280, 400]) {
      const ctx = await browser.newContext({ viewport: { width, height: 900 }, colorScheme: scheme, reducedMotion: 'reduce' });
      const p = await ctx.newPage();
      const errs = [];
      p.on('pageerror', e => errs.push(e.message));
      p.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
      await p.goto('file://' + tmp, { waitUntil: 'load' });
      await p.waitForTimeout(800);
      const sw = await p.evaluate(() => [document.documentElement.scrollWidth, window.innerWidth]);
      const file = path.join(out, `${scheme}-${width}.png`);
      await p.screenshot({ path: file, fullPage: true });
      console.log(file, 'scrollWidth/innerWidth', sw.join('/'), errs.length ? 'ERRORS: ' + errs.join(' | ') : '');
      await ctx.close();
    }
  }
  await browser.close();
  fs.unlinkSync(tmp);
})();
