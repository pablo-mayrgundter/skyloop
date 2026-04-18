import {expect, test} from '@playwright/test'

const moduleButtons = [
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

test.describe('Skyloop SPA', () => {
  test('renders main simulation and switches to module lab', async ({ page }) => {
    await page.goto('/')
    await expect(page.getByRole('heading', { name: 'Skyloop' })).toBeVisible()
    await expect(page.getByText('Simulation HUD')).toBeVisible()

    await page.getByRole('button', { name: 'Module Lab' }).click()
    await expect(page.getByText('Sections Lab')).toBeVisible()
  })

  moduleButtons.forEach((label) => {
    test(`opens module ${label}`, async ({ page }) => {
      await page.goto('/')
      await page.getByRole('button', { name: 'Module Lab' }).click()
      await page.getByRole('button', { name: label }).click()
      await expect(page.getByText(label)).toBeVisible()
    })
  })
})
