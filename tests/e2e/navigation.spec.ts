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
  for (const module of modules) {
    await page.getByRole('button', { name: module }).click()
    await expect(page.getByText(module)).toBeVisible()
  }
})
