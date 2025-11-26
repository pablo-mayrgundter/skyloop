import { test, expect } from '@playwright/test'

const modules = [
  'Sections',
  'Tethers',
  'Tube',
  'Gliders',
  'Atmosphere',
  'Weather',
  'Day/Night',
  'Bill of Materials',
  'Economics'
]

test('loads main simulation and module lab', async ({ page }) => {
  await page.goto('/')
  await expect(page.getByText('Skyloop')).toBeVisible()
  await expect(page.getByText('Main Simulation')).toBeVisible()
  await page.getByText('Module Lab').click()
  await expect(page.getByText('Modules')).toBeVisible()
})

test('module lab renders each module', async ({ page }) => {
  await page.goto('/')
  await page.getByText('Module Lab').click()
  
  // Check each module by verifying unique content appears after clicking
  const moduleChecks: Record<string, string> = {
    'Sections': 'buoyant vacuum segments',
    'Tethers': 'tether tension',
    'Tube': 'tube geometry',
    'Gliders': 'Stubbed glider physics',
    'Atmosphere': 'Atmospheric layers',
    'Weather': 'Wind profiles',
    'Day/Night': 'Day/Night Cycle',
    'Bill of Materials': 'Bill of Materials',
    'Economics': 'Economics'
  }
  
  for (const module of modules) {
    const button = page.getByRole('button', { name: module })
    await button.click()
    // Wait for the module content to appear
    await page.waitForTimeout(200)
    // Check for unique content text that appears in the module content area
    // Use locator with filter to find text containing the unique string
    const uniqueText = moduleChecks[module]
    const contentLocator = page.locator('text=' + uniqueText).first()
    await expect(contentLocator).toBeVisible({ timeout: 5000 })
  }
})
