const { chromium } = require('/opt/node22/lib/node_modules/playwright');
(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 2 });
  await page.goto('file:///home/user/Papulacandin/research_goal.html');
  const el = await page.$('#canvas');
  await el.screenshot({ path: '/home/user/Papulacandin/research_goal.png' });
  await browser.close();
  console.log('done');
})();
